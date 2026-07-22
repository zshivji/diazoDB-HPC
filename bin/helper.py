import pandas as pd
import numpy as np
import glob

def load_nif(update_index, archaea_only = False):
    if archaea_only:
        nif = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
    else:
        nif_archaea = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
        nif_bacteria = pd.read_csv('../results/hmmsearch_results/bacteria_hits.csv')
        nif = pd.concat([nif_archaea, nif_bacteria], ignore_index=True)

    nif.reset_index(inplace = True)
    nif.set_index(update_index, inplace = True) # set index to selected columns

    return nif

def filter_groups_by_unique_counts(df, group_cols, requirements):
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

    def valid_group(x):
        return all(
            x[col].nunique() >= threshold
            for col, threshold in requirements.items()
        )

    return df.groupby(group_cols).filter(valid_group)

def str_to_bool(value):
    """Convert an argument string into a boolean, supporting numbers and text."""
    # Convert to lowercase string to check string variations safely
    val_lower = str(value).lower().strip()
    
    if val_lower in ('yes', 'true', 't', 'y', '1'):
        return True
    elif val_lower in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError(f"Boolean value expected. Got '{value}'.")
