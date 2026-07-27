import pandas as pd
import argparse
import glob
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def default_proteins_dir() -> str:
    """Resolve the GTDB protein reps directory with sensible fallbacks."""
    if os.getenv("DIAZODB_PROTEIN_REPS_DIR"):
        return os.environ["DIAZODB_PROTEIN_REPS_DIR"]
    if os.path.isdir("../protein_faa_reps_latest"):
        return "../protein_faa_reps_latest"
    return "../protein_faa_reps_232"

def load_nif(update_index, archaea_only = False):
    if archaea_only:
        nif = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
    else:
        nif_archaea = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
        nif_bacteria = pd.read_csv('../results/hmmsearch_results/bacteria_hits.csv')
        nif = pd.concat([nif_archaea, nif_bacteria], ignore_index=True)

    # nif.reset_index(inplace = True)
    nif.set_index(update_index, inplace = True) # set index to selected columns

    return nif

def get_seq(genome, hit, id=None, description=None, start = None, end = None, dir: str = "../protein_faa_reps_232"):
    # get the sequence of a hit from the protein fasta file
    record_id = genome if id is None else id
    record_description = hit if description is None else description 
    file = glob.glob(f"{dir}/*/{genome}_protein.faa")[0]
    for result in SeqIO.parse(file, "fasta"):
        if result.id == hit:
            # store seq
            seq = result.seq[start:end] # trim seq to match domain len
            # convert to seqrecord
            return SeqRecord(seq, id=record_id, description=record_description)
    

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

# def str_to_bool(value):
#     """Convert an argument string into a boolean, supporting numbers and text."""
#     # Convert to lowercase string to check string variations safely
#     val_lower = str(value).lower().strip()
    
#     if val_lower in ('yes', 'true', 't', 'y', '1'):
#         return True
#     elif val_lower in ('no', 'false', 'f', 'n', '0'):
#         return False
#     else:
#         raise argparse.ArgumentTypeError(f"Boolean value expected. Got '{value}'.")
