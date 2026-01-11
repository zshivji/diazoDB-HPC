import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from itertools import combinations
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
import sys
from cluster_pos import cluster_pos

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
        
# filter for genomes to keep
genomes_to_keep.to_feather(f'../results/{dir}/nif.feather')
genomes_to_keep.to_csv(f'../results/{dir}/nif.csv')
