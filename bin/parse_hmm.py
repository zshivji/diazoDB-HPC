#!/usr/bin/env python3
"""Parse hmmsearch outputs and filter tophits results.

Inputs:
    --hits             Glob or path to hmmsearch.out results from the HMM search step.
    --outdir           Output directory for hits.feather and hits.csv.

Examples:
    python bin/parse_hmm.py \
        --hits "../results/hmmsearch_results/archaea/*.out" \
        --outdir ../results/hmmsearch_results/hits_archaea.
"""

import argparse
import glob
import re
import warnings
from pathlib import Path

import pandas as pd
from Bio import SearchIO

from cluster_pos import cluster_pos

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
    gtdb_taxonomy = pd.read_csv(
        "GTDB_taxonomy.gz",
        header=None,
        sep="\t",
        names=["GenomeID", "GTDB"],
    )

    # Create empty lists
    result_target = []
    query_id = []
    hit_id = []
    evalue = []
    best_domain_evalue = []
    bitscore = []
    bias = []
    location = []
    alens = []
    slength = []
    flag1 = []
    flag2 = []

    def append_hit(genome_id: str, gene: str, item) -> None:
        result_target.append(genome_id)
        query_id.append(gene)
        hit_id.append(item.id)
        evalue.append(item.evalue)
        best_domain_evalue.append(item.hsps[0].evalue)
        bitscore.append(item.bitscore)
        bias.append(item.bias)
        s = r"# ([0-9]+) # ([0-9]+)"
        location.append(
            re.match(s, item.description).group(1)
            + "-"
            + re.match(s, item.description).group(2)
        )
        # grab full alignment length (need to sum all domains)
        alen = 0
        for domain in item.hsps:
            alen += domain.aln_span
        alens.append(alen)
        slength.append(
            int(re.match(s, item.description).group(2))
            - int(re.match(s, item.description).group(1))
        )

    # Parse through files in output directory
    for file in glob.glob(str(hits_path)):
        # RegEx for the GenomeID (double checking that file is really a genome)
        try:
            s = r"([\w]+_[\w]+_[\d]+\.[\d])"
            genome_id = re.search(s, file).group()
        except Exception:
            continue

        # Parse file using SearchIO/HmmerIO
        for result in SearchIO.parse(file, "hmmer3-text"):
            for item in result.hits:
                # grab gene name
                s = r"([a-zA-Z]+)"  # ex. nifHDK, pchlide
                gene = re.findall(s, result.id)[0]

                # Check for positive bitscore and append the data to the corresponding lists
                if item.bitscore > 0 and item.evalue < 0.01:
                    append_hit(genome_id, gene, item)

                    # check if full seq and best domain e-val are significant
                    if item.hsps[0].evalue < 0.01:
                        flag1.append(0)
                    else:
                        # check if "full sequence Eval is sig but best domain is not"
                        flag1.append(1)

                    # check if bitscore >> bias (same order of magnitude) as bitscore
                    if item.bias != 0 and item.bitscore / item.bias > 10:
                        flag2.append(0)
                    elif item.bias == 0:
                        flag2.append(0)
                    else:
                        flag2.append(1)

    # create and store dataframe
    hits = pd.DataFrame(
        {
            "GenomeID": result_target,
            "Gene": query_id,
            "Hit": hit_id,
            "E-value": evalue,
            "Best Domain E-value": best_domain_evalue,
            "Bit Score": bitscore,
            "Bias": bias,
            "Location": location,
            "Alignment Length": alens,
            "Sequence Length": slength,
            "Flag_Eval": flag1,
            "Flag_Bias": flag2,
        }
    )

    # add taxonomy info
    hits = pd.merge(hits, gtdb_taxonomy, on="GenomeID", how="left")
    return hits


def parse_tophits(hits: pd.DataFrame, outdir: Path, min_genes: int, gene_range: int) -> None:
    # save "contig" as col
    hits["contig"] = hits["Hit"].str.split("_").str[:-1].str.join("_")

    # multi-index to cluster by genome, contig
    hits.set_index(["GenomeID", "contig"], inplace=True)
    hits.sort_index(inplace=True)
    hits.drop_duplicates(inplace=True)

    # filter for genome, contig with at least 3 unique genes (nifHDKENB)
    filtered_df = hits.groupby(level=["GenomeID", "contig"]).filter(
        lambda x: x["Gene"].nunique() >= min_genes
    )

    # make sure these 3 unique genes are not the same hit (i.e. not the same gene in reference genome)
    filtered_df2 = filtered_df.groupby(level=["GenomeID", "contig"]).filter(
        lambda x: x["Hit"].nunique() >= min_genes
    )

    genomes_to_keep = pd.DataFrame(columns=filtered_df2.columns)
    # iterate through each genome and contig
    for genome in filtered_df2.index.get_level_values(0).unique():
        for contig in filtered_df2.loc[genome].index.get_level_values(0).unique():
            tmp = filtered_df2.loc[(genome, contig)]

            # only keep numbers that have clusters >= min_genes
            pos_clusters = cluster_pos(tmp.Hit.unique(), gene_range)

            # for each cluster, find the best combination of genes (min e-value)
            for cl in pos_clusters:
                pos = [contig + "_" + str(p) for p in cl]
                no_pos = len(pos)

                # need at least min_genes genes to continue
                if no_pos < min_genes:
                    continue

                # only keep hits that are in the cluster
                tmp2 = tmp[tmp.Hit.isin(pos)].reset_index()

                genomes_to_keep = pd.concat([genomes_to_keep, tmp2])

    # export
    outdir.mkdir(parents=True, exist_ok=True)
    genomes_to_keep.to_feather(str(outdir) + "_hits.feather")
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
