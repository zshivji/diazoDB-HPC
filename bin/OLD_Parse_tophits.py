import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from itertools import combinations
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
import sys

# make sure sequences are close to each other (within 15 genes)
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
        
# make sure combos have unique genes and hits
def is_valid_combo(combo, pos):
    genes = {c[0] for c in combo}
    hits = {c[1] for c in combo}

    if len(genes) == pos and len(hits) == pos:
        return True


# get parsed hmm results
dir = sys.argv[1]
hits = pd.read_feather(f'../results/{dir}/hits.feather')

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

genomes_to_keep = []
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
            tmp2 = tmp[tmp.Hit.isin(pos)]
        
            # get all pairs of gene, hit
            items = [(gene, hit) for gene, hit in zip(tmp2.Gene, tmp2.Hit)]
    
            # get all valid combinations (max 6 possible genes--nifHDKENB)
            valid_combinations = [combo for combo in combinations(items, min(no_pos,6)) if is_valid_combo(combo, min(no_pos,6))]

            if len(valid_combinations) == 0:
                continue

            # select combo with minimum total e-value
            best_combo = []
            min_score = 1
            tmp2.set_index(['Gene', 'Hit'], append = True, inplace = True)
            
            for combo in valid_combinations:

                score = tmp2.loc[(genome, contig)].loc[list(combo), "E-value"].sum()
                if score < min_score:
                    min_score = score
                    best_combo = combo

            best_combo = [(genome, contig) + tuple(item) for item in best_combo]
            genomes_to_keep.extend(best_combo)

# filter for genomes to keep
filtered_df2.set_index(['Gene', 'Hit'], append = True, inplace = True)
filtered_df2 = filtered_df2.loc[genomes_to_keep]
filtered_df2.to_feather(f'../results/{dir}/nif.feather')
filtered_df2.to_csv(f'../results/{dir}/nif.csv')