import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

def cluster_pos(pos, limit = 15):
    # get gene position from hit
    pos_num = np.array([int(p.split('_')[-1]) for p in pos])

    # use hierarchical clustering to group neighboring genes
    Z = linkage(pos_num.reshape(-1, 1), method = 'ward')
    clusters = fcluster(Z, t = limit, criterion = 'distance')

    # only keep clusters with at least 3 genes
    cluster_no = [item for item in set(clusters) if list(clusters).count(item) >= 3]

    # store as list of lists
    filtered = []
    for cl in cluster_no:
        filtered.append(list(pos_num[clusters == cl]))
    
    return filtered