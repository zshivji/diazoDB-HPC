import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

# switch operons when 1) strand changes or 2) distance between genes is greater than limit
def cluster_operon(hits, limit = 15):

    hits = hits.copy()
    
    # sort hits by gene position and orientation
    hits["pos_num"] = hits["Hit"].str.split("_").str[-2].astype(int)
    # hits = hits.sort_values(["pos_num", "Orientation"])

    operon_counter = 1
    operon_labels = {}

    for orient, strand_hits in hits.groupby("Orientation", sort=False):
        
        # use hierarchical clustering to group neighboring genes
        strand_hits = strand_hits.sort_values("pos_num")
        pos = strand_hits["pos_num"].to_numpy().reshape(-1, 1)

        # cluster genes by distance
        Z = linkage(pos, method = 'ward') # try method = 'single'
        clusters = fcluster(Z, t = limit, criterion = 'distance')
        
        

    return clusters
