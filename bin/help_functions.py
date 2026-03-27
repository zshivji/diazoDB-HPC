import pandas as pd
import numpy as np
import glob

def load_nif(update_index, archaea_only = False):
    if archaea_only:
        nif = pd.read_feather('results/archaea/nif.feather')
    else:
        nif_archaea = pd.read_feather('results/archaea/nif.feather')
        nif_bacteria = pd.read_feather('results/bacteria/nif.feather')
        nif = pd.concat([nif_archaea, nif_bacteria], ignore_index=True)

    nif.reset_index(inplace = True)
    nif.set_index(update_index, inplace = True) # set index to selected columns

    return nif