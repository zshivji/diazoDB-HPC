import pandas as pd
from Bio import SeqIO
import glob

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from cluster_pos import cluster_pos

# make sure gene clusters have at least nifHDK
def gene_check(genes):
    if genes.__contains__('nifH'):
        if genes.__contains__('nifD'):
            if genes.__contains__('nifK'):
                return True

# multi-index to cluster by genome, contig
nif = pd.read_csv(f'nif_rescheck_nofilt.csv')
nif.reset_index(inplace = True)
nif.set_index(['GenomeID', 'contig'], inplace = True)
nif.sort_index(inplace = True)
nif.drop_duplicates(inplace = True)

# filter for genome, contig with at least 3 unique genes (nifHDKENB)
filtered_nif = nif.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Gene'].nunique() >= 3)

# make sure these 3 unique genes are not the same hit (i.e. not the same gene in reference genome)
filtered_nif2 = filtered_nif.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Hit'].nunique() >= 3)

genomes_to_keep = []
# iterate through each genome and contig
for genome in filtered_nif2.index.get_level_values(0).unique(): # iterate through each genome
    for contig in filtered_nif2.loc[genome].index.get_level_values(0).unique(): # iterate through each contig

        tmp = filtered_nif2.loc[(genome, contig)]

        # only keep numbers that have clusters >= 3
        pos_clusters = cluster_pos(tmp.Hit.unique(), 20)

        # for each cluster, find the best combination of genes (min e-value)
        for ind, cl in enumerate(pos_clusters):
            pos = [contig + '_' + str(p) for p in cl]
            no_pos = len(pos)
            
            # need at least 3 genes to continue
            if no_pos < 3:
                continue

            # only keep hits that are in the cluster
            tmp2 = tmp[tmp.Hit.isin(pos)]

            # only keep hits that are in the cluster
            tmp2 = tmp[tmp.Hit.isin(pos)].reset_index()

            # check if all genes are present
            if gene_check(tmp2.residue_match.to_list()):
                # get index
                items = [(genome, contig, hit) for hit in tmp2.Hit]
                genomes_to_keep.extend(items)

# filter for genomes to keep
filtered_nif2.set_index(['Hit'], append = True, inplace = True)
filtered_nif2 = filtered_nif2.loc[genomes_to_keep]
filtered_nif2.sort_index(inplace = True)

#clean up cols
filtered_nif2['Gene'] = filtered_nif2['residue_match']
filtered_nif2 = filtered_nif2[['Gene', 'E-value', 'Bit Score', 'Location', 'Orientation', 'Alignment Length', 'Sequence Length', 'GTDB']]
filtered_nif2.drop_duplicates(inplace = True)

# export 
filtered_nif2.to_feather(f'nif_final.feather')
filtered_nif2.to_csv(f'nif_final.csv')

# export fasta files for each gene
nif = pd.read_feather('../results/nif_final.feather')

# export fasta files (use .copy() to avoid SettingWithCopyWarning)
nifH = nif[(nif.Gene == 'nifH')].copy()
nifD = nif[(nif.Gene == 'nifD')].copy()
nifK = nif[(nif.Gene == 'nifK')].copy()
nifB = nif[(nif.Gene == 'nifB')].copy()
nifE = nif[(nif.Gene == 'nifE')].copy()
nifN = nif[(nif.Gene == 'nifN')].copy()
nifNB = nif[(nif.Gene == 'nifNB')].copy()

# anfH = nif[(nif.Gene == 'anfH')].copy()
# anfD = nif[(nif.Gene == 'anfD')].copy()
# anfK = nif[(nif.Gene == 'anfK')].copy()
# anfO = nif[(nif.Gene == 'anfO')].copy()
# anfG = nif[(nif.Gene == 'anfG')].copy()

# vnfH = nif[(nif.Gene == 'vnfH')].copy()
# vnfD = nif[(nif.Gene == 'vnfD')].copy()
# vnfK = nif[(nif.Gene == 'vnfK')].copy()
# vnfO = nif[(nif.Gene == 'vnfO')].copy()
# vnfG = nif[(nif.Gene == 'vnfG')].copy()


gene_list = [nifH, nifD, nifK, nifB, nifE, nifN, nifNB]
gene_names = ['nifH', 'nifD', 'nifK', 'nifB', 'nifE', 'nifN', 'nifNB']

# get fasta sequences for each gene & export to fasta
for gene, name in zip(gene_list, gene_names):
    print(name, flush=True)
    records = []

    for genome,hit in gene.iterrows():
        file = glob.glob(f"../all_rep_proteins_aa/*/{genome[0]}_protein.faa")[0]
        for result in SeqIO.parse(file, "fasta"):
            if result.id == genome[-1]: # -1 or 2
                # store seq
                gene.loc[genome, 'Seq'] = str(result.seq)
                # convert to seqrecord (>hit genome)
                record = SeqRecord(result.seq, id=genome[-1], description=genome[0])
                records.append(record)
                # exit loop once sequence is found
                break
    print(len(records), flush=True)    
    # Write the records to a FASTA file
    with open("../results/checked_" + name + ".fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

print("done", flush=True)