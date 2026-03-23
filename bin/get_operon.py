import pandas as pd
import os
import glob
from cluster_pos import cluster_pos

# grab checked nifHDKENB
nif = pd.read_feather('../results/nif_final_04292025.feather')

# multi-index to cluster by genome, contig
nif.reset_index(inplace = True)
nif.set_index(['GenomeID', 'contig'], inplace = True)
nif.sort_index(inplace = True)
nif.drop_duplicates(inplace = True)

# iterate through each genome and contig
for genome in nif.index.get_level_values(0).unique(): # iterate through each genome
    for contig in nif.loc[genome].index.get_level_values(0).unique(): # iterate through each contig

        tmp = nif.loc[(genome, contig)]

        # only keep numbers that have clusters >= 3
        pos_clusters = cluster_pos(tmp.Hit.unique())
        
        # for each cluster, make sure nifHDKEN (nifHDK)
        for cl in pos_clusters:
            # get positions within +/-10 genes of center
            middle = cl[len(cl)//2] # center of cluster
            pos = [middle+num for num in range(-12,13) if middle+num > 0]
            acc = [contig + '_' + str(p) for p in pos] # get acc

            # write list items into a .txt file
            with open("seqs.txt", "w") as f:
                for item in acc:
                    f.write(f"{item}\n")

            # save subsets as fasta
            file = glob.glob(f"../all_rep_proteins_aa/*/{genome}_protein.faa")[0]
            os.system(f"seqtk subseq {file} seqs.txt > ../operon-org/input-fastas/{genome}_{contig}_operon.fasta")
            os.system("rm seqs.txt")
