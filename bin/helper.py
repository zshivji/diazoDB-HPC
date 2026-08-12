import pandas as pd
import argparse
import glob
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
import Bio.Phylo.NewickIO

def default_proteins_dir() -> str:
    """Resolve the GTDB protein reps directory with sensible fallbacks."""
    if os.getenv("DIAZODB_PROTEIN_REPS_DIR"):
        return os.environ["DIAZODB_PROTEIN_REPS_DIR"]
    if os.path.isdir("../protein_faa_reps_latest"):
        return "../protein_faa_reps_latest"
    return "../protein_faa_reps_232"

def load_nif(update_index, archaea_only=False, hits_file=None):
    if hits_file:
        nif = pd.read_csv(hits_file)
    elif archaea_only:
        nif = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
    else:
        nif_archaea = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
        nif_bacteria = pd.read_csv('../results/hmmsearch_results/bacteria_hits.csv')
        nif = pd.concat([nif_archaea, nif_bacteria], ignore_index=True)

    # nif.reset_index(inplace = True)
    nif.set_index(update_index, inplace = True) # set index to selected columns

    return nif

def get_seq(genome, hit, id=None, description=None, start=None, end=None, dir: str = "../protein_faa_reps_232"):
    # get the sequence of a hit from the protein fasta file
    record_id = genome if id is None else id
    record_description = hit if description is None else description
    candidates = glob.glob(f"{dir}/*/{genome}_protein.faa")
    if not candidates:
        raise FileNotFoundError(
            f"Protein FASTA not found for genome {genome!r} under {dir}"
        )
    file = candidates[0]
    for result in SeqIO.parse(file, "fasta"):
        if result.id == hit:
            # store seq
            seq = result.seq[start:end] # trim seq to match domain len
            # convert to seqrecord
            return SeqRecord(seq, id=record_id, description=record_description)
    raise KeyError(f"Sequence {hit!r} was not found in {file}")

def filter_groups_by_unique_counts(df, group_cols, requirements, exclude='nifB'):
    """
    Filter groups based on multiple unique-count requirements.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    group_cols : list[str]
        Columns defining groups
    requirements : dict
        Dictionary of {column: minimum_unique_count}

    Example
    -------
    requirements={
        "Gene": 3,
        "Hit": 3
    }
    """

    def has_excluded_protein(x):
        return (
            exclude is not None
            and 'Gene' in x.columns
            and exclude in x['Gene'].unique()
        )

    def meets_requirements(x):
        return all(
            x[col].nunique() >= threshold
            for col, threshold in requirements.items()
        )

    # Step 1: remove genomes that fewer than 3 unique nif genes.
    df = df.copy()
    df = df.groupby('GenomeID').filter(
        lambda x: meets_requirements(x)
    )

    # Step 2: remove groups (typically operons) that do not meet the requirements or contain the excluded protein.
        # an operon can be kept with fewer than (3) unique nif genes if it contains the excluded protein (e.g., nifB)
    return df.groupby(group_cols).filter(
        lambda x: has_excluded_protein(x) or meets_requirements(x)
    )

def tree_node_match_metadata(tree_file):

    clusters = pd.read_csv('../results/final/nif_clusters.csv')

    tree = Phylo.read(tree_file, "newick")

    for clade in tree.get_terminals():
        # for each node, update nodeID to match the metadata headers (except added outbranches)
        if clade.name.startswith('WP_'):
            new_name = clade.name
        else:
            protein = clade.name
            contig = '_'.join(protein.split('_')[:-1])
            position = protein.split('_')[-1]
            cluster = clusters[clusters["contig"] == contig]
            for _, row in cluster.iterrows():
                positions = str(row['pos_num']).strip('[]').replace(',', ' ').split()
                if position in {p.strip() for p in positions}:
                    organism = row['Organism']
                    clusterID = row['cluster']
                    genome = row['GenomeID']
                    contig = row['contig']
                    operon = row['operon']
                    new_name = f"{organism} | {clusterID} | {genome} | {contig} | {operon}"
                    break

        clade.name = new_name

    Phylo.write(tree, tree_file, "newick")


def main():
    parser = argparse.ArgumentParser(description="DiazoDB helper functions")
    subparsers = parser.add_subparsers(dest="command", required=True)

    tree_parser = subparsers.add_parser(
        "tree_node_match_metadata",
        help="Replace tree terminal IDs with metadata-matched IDs.",
    )
    tree_parser.add_argument("tree_file", help="Newick tree file to update in place.")

    args = parser.parse_args()
    if args.command == "tree_node_match_metadata":
        tree_node_match_metadata(args.tree_file)


if __name__ == "__main__":
    main()
