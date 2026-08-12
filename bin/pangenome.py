#!/usr/bin/env python3
"""Parse microbannotator results and visualize with pangenome figure.

Inputs:
  

Examples:
    python bin/pangenome.py 
    
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

def parse_annotations() -> pd.DataFrame: #hits_path: Path, outdir: Path
    gtdb_taxonomy = pd.read_csv(
        "GTDB_taxonomy.gz",
        header=None,
        sep="\t",
        names=["GenomeID", "GTDB"],
    )

    annotation_results = pd.DataFrame(columns=["genome_id", "query_id", "protein_id", "product", "ko_number", "taxonomy",
                                               "function_go", "compartment_go", "process_go", "interpro",
                                               "pfam", "ec_number", "database"])
    # Parse through files in output directory
    for file in glob.glob(str("../pangenome/batch_*/annotation_results/*.annot")):
        # add annotation results to the dataframe
        tmp = pd.read_csv(file, sep="\t")
        tmp["genome_id"] = re.search(r'/(*?)_[0-9]_[a-z]{3}_operon.fasta.annot' ,file).group(0)
        print(tmp["genome_id"][0])
        annotation_results = pd.concat([annotation_results, tmp], ignore_index=True)

    annotations = list(set(annotation_results['ko_number'].to_list()))
    annotations.extend('genomeID')

    heat_map = pd.DataFrame(columns=annotations)

    annotation_results.set_index('genome_id', inplace=True)
    for genome in annotation_results.index.unique():
        genome_annotations = annotation_results.loc[genome]
        heat_map.loc[genome] = genome_annotations['ko_number'].value_counts()

    heat_map.fillna(0, inplace=True)
    heat_map.to_csv("../pangenome/heatmap.csv")



def main() -> None:
    parse_annotations()
    # args = parse_args()

    # if args.hits and args.outdir:
    #     hits = parse_hits(Path(args.hits), Path(args.outdir))
    #     parse_tophits(hits, Path(args.outdir), args.min_genes, args.gene_range)
    # elif args.hits or args.outdir:
    #     raise SystemExit("Both --hits and --outdir are required to parse HMM hits."
    #                      "Please provide both arguments and try again." )

if __name__ == "__main__":
    main()
