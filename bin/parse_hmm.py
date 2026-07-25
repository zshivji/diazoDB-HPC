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
import re
import warnings
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SearchIO
from scipy.cluster.hierarchy import linkage, fcluster
from helper import filter_groups_by_unique_counts

warnings.filterwarnings("ignore", category=FutureWarning)


def genome_id_from_path(file: Path) -> str | None:
    """Extract a GTDB accession or runner input name from a domtblout path."""
    accession_match = re.search(r"([\w]+_[\w]+_[\d]+\.[\d])", str(file))
    if accession_match:
        return accession_match.group()
    if "__job_" in file.name:
        return file.name.split("__job_", 1)[0]
    return None


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
    print(f"Parsing HMM hits...\n\t\t from {hits_path}\n\t\t saving results to {outdir}_hits.csv", flush=True)

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

    ### NEED TO UPDATE TO ACCOMODATE ALL INPUT TYPES
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
        s = r"# ([0-9]+) # ([0-9]+) # (-?1) #"  ### NEED TO UPDATE TO ACCOMODATE ALL INPUT TYPES
        # location = full protein length
        location.append(
            re.match(s, hit.description).group(1)
            + "-"
            + re.match(s, hit.description).group(2) 
        )
        orientation.append(re.match(s, hit.description).group(3))
 

    # Parse through files in output directory
    for file in hits_path.glob("*.domtblout"):
        # Prefer a GTDB accession when present. Runner intermediates use
        # <input-stem>__job_<uuid>_nif.domtblout, so uploaded files without a
        # GTDB accession fall back to their sanitized input stem.
        genome_id = genome_id_from_path(file)
        if genome_id is None:
            print(f"Could not determine GenomeID from {file}", flush=True)
            continue

        # Parse file using SearchIO/HmmerIO w/ --domtblout option
        for result in SearchIO.parse(file, "hmmsearch3-domtab"):
            # grab gene name, ex. nifHDK, pchlide
            gene = re.findall(r"([a-zA-Z]+)", result.id)[0]
            
            for hit in result.hits:
                # Check for positive bitscore and append the data to the corresponding lists
                if hit.bitscore <= 0 or hit.evalue >= 0.001:
                    continue

                for hsp in hit.hsps:
                    if hsp.bitscore <= 0 or hsp.evalue >= 0.001:
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

    # pre-filter to remove contigs strands with less than min_genes unique nif genes)
    hits = filter_groups_by_unique_counts(
        hits,
        group_cols=["GenomeID", "contig", "Orientation"],
        requirements={
        "Gene": min_genes,
        "Hit": min_genes
        },
        exclude='nifB'
    )

    # store genomes meeting criteria in new dataframe
    genomes_to_keep = pd.DataFrame(columns=hits.columns)

    # iterate through each genome and contig and split into operons
    for (genome, contig), tmp in hits.groupby(["GenomeID", "contig"]):

            # split contigs into operons when 1) strand changes or 2) distance between genes is greater than limit
            operon_counter = 1

            # sort hits by gene position and orientation
            # for each strand, cluster genes by distance and assign operon labels
            for _, strand in tmp.groupby("Orientation"):

                if len(strand) < 2: # skip strands with only nifB (linkage clustering breaks)
                    clusters = np.array([operon_counter])
                else: 
                    # use hierarchical clustering to group neighboring genes
                    strand = strand.sort_values("pos_num")
                    pos = strand["pos_num"].to_numpy().reshape(-1, 1)

                    # cluster genes by distance
                    Z = linkage(pos, method = 'single') # try method = 'ward' but singe prob better
                    clusters = fcluster(Z, t = gene_range, criterion = 'distance')
                    clusters = clusters.astype(int) + operon_counter -1 # convert to int for easier handling

                strand["operon"] = clusters

                # save clusters to tmp2
                operon_counter += clusters.max()
                genomes_to_keep = pd.concat([genomes_to_keep, strand])

    # refilter to only keep clusters with at least 3 unique genes, unless a gene is excluded from the requirement (ex. nifB)
    genomes_to_keep = filter_groups_by_unique_counts(
        genomes_to_keep,
        group_cols=["GenomeID", "contig", "operon"],
        requirements={
        "Gene": min_genes,
        "Hit": min_genes
        },
        exclude='nifB'
    )

    # record which gene has the highest bitscore for each protein
    genomes_to_keep["top_hit"] = (genomes_to_keep.sort_values("Bit Score", ascending=False).groupby(["GenomeID", "contig", "Hit"])["Gene"].transform("first"))

    # add taxonomy info (make sure metadata is the same release as the GTDB rep seqs)
    ### NNED TP UPDATE GTDB INPUT PATH
    gtdb_taxonomy = pd.read_csv("GTDB_metadata.gz", sep="\t", usecols=['accession', 'gtdb_taxonomy']).rename(columns={"accession": "GenomeID", "gtdb_taxonomy": "GTDB"})
    genomes_to_keep = pd.merge(genomes_to_keep, gtdb_taxonomy, on="GenomeID", how="left")

    # export
    outdir.mkdir(parents=True, exist_ok=True)
    genomes_to_keep.to_csv(str(outdir) + "_hits.csv", index=False)
    print(f"Saved {genomes_to_keep.shape[0]} hits to {outdir}_hits.csv", flush=True)


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
