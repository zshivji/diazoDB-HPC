"""Filter HMM hits to clustered contigs and write top hits outputs.

Inputs:
    --hits    Path to hits.feather produced by the HMM parsing step.
    --outdir  Output directory for tophits.feather and tophits.csv.

Example:
    python bin/Parse_tophits.py \
        --hits ../results/run_2026_05_22/hits.feather \
        --outdir ../results/run_2026_05_22
"""

import argparse
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
from itertools import combinations
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
import argparse
from cluster_pos import cluster_pos
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter HMM hits to top contig clusters with >= 3 unique genes. "
            "Writes tophits.feather and tophits.csv to the output directory."
        )
    )
    parser.add_argument(
        "--hits",
        required=True,
        help="Path to hits.feather from the HMM parsing step.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for tophits.feather and tophits.csv.",
    )
    parser.add_argument(
        "--reload-fasta",
        action="store_true",
        help=(
            "Reload FASTA sequences for downstream steps if needed. "
            "This flag is accepted for pipeline compatibility but is not used here."
        ),
    )
    return parser.parse_args()


def parse_tophits(hits_path: Path, outdir: Path) -> None:
    hits = pd.read_feather(hits_path)

    # save "contig" as col
    hits['contig'] = hits['Hit'].str.split('_').str[:-1].str.join('_') 

    # multi-index to cluster by genome, contig
    hits.set_index(['GenomeID', 'contig'], inplace = True)
    hits.sort_index(inplace = True)
    hits.drop_duplicates(inplace = True)

    # filter for genome, contig with at least 3 unique genes (nifHDKENB)
    filtered_df = hits.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Gene'].nunique() >= 3)

    # make sure these 3 unique genes are not the same hit (i.e. not the same gene in reference genome)
    filtered_df2 = filtered_df.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Hit'].nunique() >= 3)

    genomes_to_keep = pd.DataFrame(columns = filtered_df2.columns)
    # iterate through each genome and contig
    for genome in filtered_df2.index.get_level_values(0).unique(): # iterate through each genome
        for contig in filtered_df2.loc[genome].index.get_level_values(0).unique(): # iterate through each contig

            tmp = filtered_df2.loc[(genome, contig)]

            # only keep numbers that have clusters >= 3
            pos_clusters = cluster_pos(tmp.Hit.unique())

            # for each cluster, find the best combination of genes (min e-value)
            for cl in pos_clusters:
                pos = [contig + '_' + str(p) for p in cl]
                no_pos = len(pos)
                
                # need at least 3 genes to continue
                if no_pos < 3:
                    continue

                # only keep hits that are in the cluster
                tmp2 = tmp[tmp.Hit.isin(pos)].reset_index()

                genomes_to_keep = pd.concat([genomes_to_keep, tmp2])
    
    # expot
    outdir.mkdir(parents=True, exist_ok=True)
    genomes_to_keep.to_feather(outdir / "tophits.feather")
    genomes_to_keep.to_csv(outdir / "tophits.csv", index=False)


def main() -> None:
    warnings.filterwarnings("ignore", category=FutureWarning)
    args = parse_args()
    parse_tophits(Path(args.hits), Path(args.outdir))


if __name__ == "__main__":
    main()