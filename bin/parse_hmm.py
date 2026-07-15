#!/usr/bin/env python3
"""Parse hmmsearch outputs and filter tophits results.

Inputs:
    --hits             Glob or path to hmmsearch.out results from the HMM search step.
    --outdir           Output directory for hits.feather and hits.csv.

Examples:
    python bin/parse_hmm.py \
        --hits ../results/hmmsearch_results/archaea/ \
        --outdir ../results/hmmsearch_results/archaea
"""

import argparse
import glob
import re
import warnings
from pathlib import Path
import os

import pandas as pd
import numpy as np
from Bio import SearchIO

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from helper import filter_groups_by_unique_counts

warnings.filterwarnings("ignore", category=FutureWarning)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse hmmsearch outputs and filter tophits. "
            "Writes hits.feather/hits.csv and tophits.feather/tophits.csv."
        )
    )
    parser.add_argument(
        "--hits",
        help="Path or glob to hmmsearch.out results from the HMM search step.",
    )
    parser.add_argument(
        "--outdir",
        help="Output directory for hmm hits that meet the following criteria (1) significant evalue, (2) have at least x nif genes within a y gene range.",
    )

    parser.add_argument(
        "--min_genes",
        type=int,
        default=3,
        help="Minimum number of nif genes required to be considered a tophit cluster (default: 3)."
    )
    
    parser.add_argument(
        "--gene_range",
        type=int,
        default=15,
        help="Maximum gene range (+/- from gene of interest) for nif genes to be considered a cluster (default: 15).")
    
    return parser.parse_args()


# Get taxonomy from GTDB, https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/
# NCBI taxonomy from metadata files

def parse_hits(hits_path: Path, outdir: Path) -> pd.DataFrame:
    print(f"Parsing HMM hits...\n\t\t from {hits_path}\n\t\t saving results to {outdir}_hits.csv")

    # Create empty lists
    result_target = []
    query_id = []
    hit_id = []
    evalue = []
    bitscore = []
    bias = []
    dom_start = []
    dom_end = []
    location = []
    orientation = []
    flag = []

    ind_map = 'ABCDEFGHIJKLM'
    def append_hit(genome_id: str, gene: str, hit, hsp) -> None:
        # Save hit domains and relevant metadata
        result_target.append(genome_id)
        query_id.append(gene)
        hit_id.append(hit.id + "_" + ind_map[hsp.domain_index-1]) # hit in the form of contig_geneNo_domA/B/C
        evalue.append(hsp.evalue) # domain independent evalue
        bitscore.append(hsp.bitscore) # domain bitscore
        bias.append(hsp.bias) # domain bias
        dom_start.append(hsp.env_start) 
        dom_end.append(hsp.env_end)

        # grab description parameters
        s = r"# ([0-9]+) # ([0-9]+) # (-?1) #"
        # location = full protein length
        location.append(
            re.match(s, hit.description).group(1)
            + "-"
            + re.match(s, hit.description).group(2) 
        )
        orientation.append(re.match(s, hit.description).group(3))
 

    # Parse through files in output directory
    for file in hits_path.glob("*.domtblout"):
        # RegEx for the GenomeID (double checking that file is really a genome)
        try:
            s = r"([\w]+_[\w]+_[\d]+\.[\d])"
            genome_id = re.search(s, str(file)).group()
        except Exception as e:
            print(f"An error occurred: {e}")
            continue
            
        # Parse file using SearchIO/HmmerIO w/ --domtblout option
        for result in SearchIO.parse(file, "hmmsearch3-domtab"):
            # grab gene name, ex. nifHDK, pchlide
            gene = re.findall(r"([a-zA-Z]+)", result.id)[0]
            
            for hit in result.hits:
                # Check for positive bitscore and append the data to the corresponding lists
                if hit.bitscore <= 0 and hit.evalue >= 0.01:
                    continue

                for hsp in hit.hsps:
                    if hsp.bitscore <= 0 and hsp.evalue >= 0.01:
                        continue
                    append_hit(genome_id, gene, hit, hsp)

                    # check if bitscore >> bias (same order of magnitude) as bitscore
                    if hit.bias != 0 and hit.bitscore / hit.bias > 10:
                        flag.append(0)
                    elif hit.bias == 0:
                        flag.append(0)
                    else:
                        flag.append(1)

    # create and store dataframe
    hits = pd.DataFrame(
        {
            "GenomeID": result_target,
            "Gene": query_id,
            "Hit": hit_id,
            "E-value": evalue,
            "Bit Score": bitscore,
            "Bias": bias,
            "Domain Start": dom_start,
            "Domain End": dom_end,
            "Location": location,
            "Orientation": orientation,
            "Flag_Bias": flag
        }
    )

    hits.drop_duplicates(inplace=True)
    return hits

def parse_tophits(hits: pd.DataFrame, outdir: Path, min_genes: int, gene_range: int) -> None:
    # save "contig" and "pos_num" as col
    hits["contig"] = hits["Hit"].str.split("_").str[:-2].str.join("_")
    hits['pos_num'] = hits["Hit"].str.split("_").str[-2].astype(int)
    hits['operon'] = np.nan

    # pre-filter to remove contigs strands with less than min_genes unique nif genes
    hits = filter_groups_by_unique_counts(
        hits,
        group_cols=["GenomeID", "contig", "Orientation"],
        requirements={
        "Gene": 3,
        "Hit": 3
        }
    )

    # store genomes meeting criteria in new dataframe
    genomes_to_keep = pd.DataFrame(columns=hits.columns)

    # iterate through each genome and contig to check criteria
        # each operon has at least 3 unique nif genes within a gene range of 15 (default)
    for (genome, contig), tmp in hits.groupby(["GenomeID", "contig"]):

            # split contigs into operons when 1) strand changes or 2) distance between genes is greater than limit
            operon_counter = 1
            
            # sort hits by gene position and orientation
            # for each strand, cluster genes by distance and assign operon labels
            for _, strand in tmp.groupby("Orientation"):
                
                # use hierarchical clustering to group neighboring genes
                strand = strand.sort_values("pos_num")
                pos = strand["pos_num"].to_numpy().reshape(-1, 1)

                # cluster genes by distance
                Z = linkage(pos, method = 'single') # try method = 'ward' but singe prob better
                clusters = fcluster(Z, t = gene_range, criterion = 'distance')
                clusters = clusters.astype(int) + operon_counter -1 # convert to int for easier handling

                strand["operon"] = clusters

                # for each cluster, only keep clusters with at least 3 unique genes
                strand = filter_groups_by_unique_counts(
                    strand,
                    group_cols=["GenomeID", "contig", "operon"],
                    requirements={
                    "Gene": 3,
                    "Hit": 3
                    }
                )

                # save clusters to tmp2
                operon_counter += clusters.max()
                genomes_to_keep = pd.concat([genomes_to_keep, strand])

    # add taxonomy info
    gtdb_taxonomy = pd.read_csv(
        "GTDB_taxonomy.gz",
        header=None,
        sep="\t",
        names=["GenomeID", "GTDB"],
    )
    
    genomes_to_keep = pd.merge(genomes_to_keep, gtdb_taxonomy, on="GenomeID", how="left")

    # export
    outdir.mkdir(parents=True, exist_ok=True)
    genomes_to_keep.to_csv(str(outdir) + "_hits.csv", index=False)


def main() -> None:
    args = parse_args()

    if args.hits and args.outdir:
        hits = parse_hits(Path(args.hits), Path(args.outdir))
        parse_tophits(hits, Path(args.outdir), args.min_genes, args.gene_range)
    elif args.hits or args.outdir:
        raise SystemExit("Both --hits and --outdir are required to parse HMM hits."
                         "Please provide both arguments and try again." )

if __name__ == "__main__":
    main()
