from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import numpy as np
import os
import glob
import sys

date = sys.argv[1]

def check_gene(gene, ref_seq, important_residues, passing_score, p=False):
    
    alignment = AlignIO.read(f"../results/{gene}.aln", "fasta")

    for record in alignment:
        if ref_seq in record.description:
            aln = record.seq
            break

    # map important cols in aln to ref_seq
    residue_to_alignment = {}
    residue_idx = 0  # index in original (ungapped) sequence
    col2res = {}

    for aln_idx, char in enumerate(aln):
        if char != '-':
            residue_idx += 1
            if residue_idx in important_residues:
                residue_to_alignment[residue_idx] = aln_idx
                col2res[aln_idx] = aln[aln_idx]
                
            # early stop if we've found everything
            if len(residue_to_alignment) == len(important_residues):
                break

    # print the corresponding residues in the original sequence
    if p:
        for residue, aln_idx in residue_to_alignment.items():
            print(f"Residue {residue} corresponds to alignment index {aln_idx}: {aln[aln_idx]}", flush=True)

    # store in dataframe        
    acc = [result.id for result in alignment]
    seqs = [list(str(result.seq)) for result in alignment]
    hits = [result.description.split(' ')[-1] for result in alignment]

    pssm = pd.DataFrame(seqs, index = acc)

    pssm = pssm.iloc[:, list(residue_to_alignment.values())]
    pssm['hit'] = hits
    pssm['contig'] = pssm['hit'].str.split('_').str[:-1].str.join('_')

    # check if cols contain correct residue for function
    def check_res(row):
        score = 0
        for col in residue_to_alignment.values():
            if row[col] == aln[col]:
                score += 2
                if aln[col] == 'C': # higher weight for correct C
                    score += 1
            elif aln[col] == 'C': # if ref seq is C, must also be C
                continue
            elif row[col] != '-': # greater penalty for gap than incorrect residue
                score += 1
        if score >= passing_score:
            return score
        else:
            return np.nan
        #return score

    pssm['score'] = pssm.apply(check_res, axis = 1)
    pssm.dropna(subset = ['score'], inplace = True)
    pssm.rename(columns=col2res, inplace=True)

    return pssm

print("Checking for conserved residues...", flush=True)
nif_archaea = pd.read_feather('../results/archaea/nif.feather')
nif_bacteria = pd.read_feather('../results/bacteria/nif.feather')
nif = pd.concat([nif_archaea, nif_bacteria])

nif.reset_index(inplace = True)
nif.set_index(['GenomeID'], inplace = True)
nif['Seq'] = ''

# nifH
print("checking nifH", flush=True)
# grab aln numbering to map to ref seq (A. vinelandii)
important_residues_nifH = [10, 11, 12, 13, 15, 16, 17, 40, 42, 44, 98, 130, 133] # azotobacter
#GKGGGKS DKD CDC --> also returns BchL, ChlL, and BchX

ref_seq_nifH = 'NZ_BSFG01000004.1_46'

passing_score = 25

nifH_checked = check_gene('nifH_split.00001', ref_seq_nifH, important_residues_nifH, passing_score, p=True)

for file in glob.glob(f'../results/nifH_split.00*.aln'):
    if file == '../results/nifH_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifH, important_residues_nifH, passing_score)
    nifH_checked = pd.concat([nifH_checked, new])

    nifH_checked.drop_duplicates(subset = ['hit'], inplace = True)
    nifH_checked.set_index(['hit'], append= True, inplace = True)

print(str(len(nifH_checked)) + " nifH seqs")

# nifD
print('checking nifD', flush=True)
important_residues_nifD = [62, 88, 154, 191, 195, 275, 442] # azotobacter
ref_seq_nifD = 'NZ_BSFG01000004.1_45' # azo
passing_score = 15

nifD_checked = check_gene('nifD_split.00001', ref_seq_nifD, important_residues_nifD, passing_score, p=True)

for file in glob.glob(f'../results/nifD_split.00*.aln'):
    if file == '../results/nifD_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifD, important_residues_nifD, passing_score)
    nifD_checked = pd.concat([nifD_checked, new])

# for each genome, only keep the best hit per gene cluster
nifD_checked.drop_duplicates(subset = ['hit'], inplace = True)
nifD_checked.set_index(['hit'], append= True, inplace = True)
nifD_checked['gene_cluster'] = nif.loc[(nifD_checked.index.get_level_values(0), nifD_checked.index.get_level_values(1), 'nifD'), 
                                       'gene_cluster'].values
nifD_checked = nifD_checked.loc[nifD_checked.groupby(['contig', 'gene_cluster'])['score'].idxmax()]

print(str(len(nifD_checked)) + " nifD seqs", flush=True)

# nifE 
print('checking nifE', flush=True)
important_residues_nifE = [25, 37, 62, 124] # azo, known catalytic
ref_seq_nifE = 'NZ_BSFG01000004.1_39' # azo

passing_score = 12 # could be more rigorous

nifE_checked = check_gene('nifE_split.00001', ref_seq_nifE, important_residues_nifE, passing_score, p=True)

for file in glob.glob(f'../results/nifE_split.00*.aln'):
    if file == '../results/nifE_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifE, important_residues_nifE, passing_score)
    nifE_checked = pd.concat([nifE_checked, new])

nifE_checked.drop_duplicates(subset = ['hit'], inplace = True)
nifE_checked.set_index(['hit'], append= True, inplace = True)

#remove nifD hits from nifE
nifE_checked = nifE_checked[~nifE_checked.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())]

#for each genome, only keep the best hit per gene cluster
nifE_checked['gene_cluster'] = nif.loc[(nifE_checked.index.get_level_values(0), nifE_checked.index.get_level_values(1), 'nifE'), 
                                       'gene_cluster'].values
nifE_checked = nifE_checked.loc[nifE_checked.groupby(['contig', 'gene_cluster'])['score'].idxmax()]

print(str(len(nifE_checked)) + " nifE seqs", flush=True)

# nifK
print('checking nifK', flush=True)
important_residues_nifK = [70, 95, 153, 188] # azotobacter
ref_seq_nifK = 'NZ_BSFG01000004.1_44'
passing_score = 9 # S188 has not been verified to be conserved (8 is too conservative)

nifK_checked = check_gene('nifK_split.00001', ref_seq_nifK, important_residues_nifK, passing_score, p=True)

for file in glob.glob(f'../results/nifK_split.00*.aln'):
    if file == '../results/nifK_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifK, important_residues_nifK, passing_score)
    nifK_checked = pd.concat([nifK_checked, new])

nifK_checked.drop_duplicates(subset = ['hit'], inplace = True)
nifK_checked.set_index(['hit'], append= True, inplace = True)

#remove nifD,E hits from nifK
nifK_checked = nifK_checked[~nifK_checked.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())]
nifK_checked = nifK_checked[~nifK_checked.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())]

#for each genome, only keep the best hit per gene cluster
nifK_checked['gene_cluster'] = nif.loc[(nifK_checked.index.get_level_values(0), nifK_checked.index.get_level_values(1), 'nifK'), 
                                       'gene_cluster'].values
nifK_checked = nifK_checked.loc[nifK_checked.groupby(['contig', 'gene_cluster'])['score'].idxmax()]

print(str(len(nifK_checked)) + " nifK seqs", flush=True)

# nifN
print('checking nifN', flush=True)
important_residues_nifN = [44] # azotobacter
ref_seq_nifN = 'NZ_BSFG01000004.1_38'

passing_score = 3

nifN_checked = check_gene('nifN_split.00001', ref_seq_nifN, important_residues_nifN, passing_score, p=True)

for file in glob.glob(f'../results/nifN_split.00*.aln'):
    if file == '../results/nifN_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifN, important_residues_nifN, passing_score)
    nifN_checked = pd.concat([nifN_checked, new])

nifN_checked.drop_duplicates(subset = ['hit'], inplace = True)
nifN_checked.set_index(['hit'], append= True, inplace = True)

# remove nifD,E,K hits from nifN
nifN_checked = nifN_checked[~nifN_checked.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())]
nifN_checked = nifN_checked[~nifN_checked.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())]
nifN_checked = nifN_checked[~nifN_checked.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())]

# for each genome, only keep the best hit per gene cluster
nifN_checked['gene_cluster'] = nif.loc[(nifN_checked.index.get_level_values(0), nifN_checked.index.get_level_values(1), 'nifN'), 
                                        'gene_cluster'].values
nifN_checked = nifN_checked.loc[nifN_checked.groupby(['contig', 'gene_cluster'])['score'].idxmax()]


print(str(len(nifN_checked)) + " nifN seqs", flush=True)

#get seq that failed checks
print('getting failed sequences', flush=True)

seqs = []

for gene in 'DKEN':

    # for each gene, get all hmm hit acc
    result = list(SeqIO.parse(f"../results/nif{gene}.fasta", "fasta"))
    hit = [record.description.split(" ")[-1] for record in result]

    # get seq that failed check
    for record, acc in zip(result, hit):
        if acc not in eval(f"list(nif{gene}_checked.index.get_level_values(1).unique())"):
            seq = SeqRecord(Seq(record.seq), id=record.id, description=acc)
            seqs.append(record)

print(str(len(seqs)) + " seqs failed nifDKEN checks", flush=True)

# Write the records to a FASTA file
with open(f"../results/nifDKEN_{date}.fasta", "w") as output_handle:
    SeqIO.write(seqs, output_handle, "fasta")

# add reference sequences
for gene in 'DKEN':
    os.system(f"seqtk subseq ../results/nif{gene}.fasta ../results/ref_seq.ids >> ../results/nifDKEN_{date}.fasta")

print('aligning failed sequences', flush=True)
# aln all seqs
num = int(len(seqs)/200) +1  # how many splits
os.system(f"seqtk split -n {num} ../results/fasta_splits/nifDKEN_split ../results/nifDKEN_{date}.fasta") # split fasta file
for i in range(num):
    print(i+1)
    os.system(f"seqtk subseq ../results/nifDKEN_{date}.fasta ../results/ref_seq.ids >> ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.fa") # add reference sequences
    os.system(f"mafft --auto --quiet --thread 4 ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.fa > ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.aln")
    

# backup check
print('checking failed sequences', flush=True)

nifD_backup = check_gene('nifDKEN_split.00001', ref_seq_nifD, important_residues_nifD, passing_score)
nifK_backup = check_gene('nifDKEN_split.00001', ref_seq_nifK, important_residues_nifK, passing_score)
nifE_backup = check_gene('nifDKEN_split.00001', ref_seq_nifE, important_residues_nifE, passing_score)
nifN_backup = check_gene('nifDKEN_split.00001', ref_seq_nifN, important_residues_nifN, passing_score)

for file in glob.glob(f'nifDKEN_split.00*.aln'):
    if file == '../results/nifDKEN_split.00001.aln':
        continue
    nifD_backup = check_gene(file[:-4], ref_seq_nifD, important_residues_nifD, passing_score)
    nifD_backup = pd.concat([nifD_backup, new])

    nifK_backup = check_gene(file[:-4], ref_seq_nifK, important_residues_nifK, passing_score)
    nifK_backup = pd.concat([nifK_backup, new])

    nifE_backup = check_gene(file[:-4], ref_seq_nifE, important_residues_nifE, passing_score)
    nifE_backup = pd.concat([nifE_backup, new])

    nifN_backup = check_gene(file[:-4], ref_seq_nifN, important_residues_nifN, passing_score)
    nifN_backup = pd.concat([nifN_backup, new])

# set index as genome, hit
nifD_backup.set_index(['hit'], append= True, inplace = True)

# remove hits that are already in saved
nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# for each genome, only keep the best hit per contig
nifD_backup['gene_cluster'] = 0
nifD_backup = nifD_backup.loc[nifD_backup.groupby(['contig'])['score'].idxmax()]
nifD_backup.drop_duplicates(inplace = True)

print(str(len(nifD_backup.index.unique())) + " nifD seqs", flush=True)
nifD_checked = pd.concat([nifD_checked, nifD_backup])

# set index as genome, hit
nifK_backup.set_index(['hit'], append= True, inplace = True)

# remove hits that are already in saved
nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# for each genome, only keep the best hit per contig
nifK_backup['gene_cluster'] = 0
nifK_backup = nifK_backup.loc[nifK_backup.groupby(['contig'])['score'].idxmax()]
nifK_backup.drop_duplicates(inplace = True)

print(str(len(nifK_backup.index.unique())) + " nifK seqs", flush=True)
nifK_checked = pd.concat([nifK_checked, nifK_backup])

# set index as genome, hit
nifE_backup.set_index(['hit'], append= True, inplace = True)

# remove hits that are already in saved
nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# for each genome, only keep the best hit per gene cluster
nifE_backup['gene_cluster'] = 0
nifE_backup = nifE_backup.loc[nifE_backup.groupby(['contig'])['score'].idxmax()]
nifE_backup.drop_duplicates(inplace = True)

print(str(len(nifE_backup.index.unique())) + " nifE seqs", flush=True)
nifE_checked = pd.concat([nifE_checked, nifE_backup])

# set index as genome, hit
nifN_backup.set_index(['hit'], append= True, inplace = True)

# remove hits that are already in saved
nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# for each genome, only keep the best hit per gene cluster
nifN_backup['gene_cluster'] = 0
nifN_backup = nifN_backup.loc[nifN_backup.groupby(['contig'])['score'].idxmax()]
nifN_backup.drop_duplicates(inplace = True)

print(str(len(nifN_backup.index.unique())) + " nifN seqs", flush=True)
nifN_checked = pd.concat([nifN_checked, nifN_backup])

# append updated annotation (based on conserved residue matching) to nif files
nif['residue_match'] = ''
nif['backup_match'] = ''

nif.set_index(['Hit'], append = True, inplace = True)
nif.sort_index(inplace = True)

# update residue match column in nif df
for gene in 'HDKEN':
    for genome, cols in eval(f"nif{gene}_checked.iterrows()"):
        nif.loc[(genome, cols.hit), 'residue_match'] = "nif" + gene

 # add backup check       
for gene in 'DKEN':
    for genome, cols in eval(f"nif{gene}_backup.iterrows()"):
        nif.loc[(genome[0], genome[1]), 'backup_match'] = "nif" + gene

# filter to get hits that passed residue matching
nifH = nif[(nif.residue_match == 'nifH') & (nif.Gene == 'nifH') & (nif['Alignment Length'] > 200)]
# nifB = nif[(nif.Gene == 'nifB')] # not done
# nifB['residue_match'] = 'nifB'

# only index OG matches
nifD = nif[((nif.residue_match == 'nifD') & (nif.Gene == 'nifD') & (nif.backup_match != 'nifD') & (nif['Alignment Length'] > 300))]
# add backup check
nifD_backup = nif[(nif.backup_match == 'nifD') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
nifD_backup = nifD_backup.groupby(['GenomeID', 'Hit']).first()
nifD = pd.concat([nifD, nifD_backup])

# only index OG matches
nifK = nif[((nif.residue_match == 'nifK') & (nif.Gene == 'nifK') & (nif.backup_match != 'nifK') & (nif['Alignment Length'] > 300))]
# add backup check
nifK_backup = nif[(nif.backup_match == 'nifK') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
nifK_backup = nifK_backup.groupby(['GenomeID', 'Hit']).first()
nifK = pd.concat([nifK, nifK_backup])

# only index OG matches
nifE = nif[((nif.residue_match == 'nifE') & (nif.Gene == 'nifE') & (nif.backup_match != 'nifE') & (nif['Alignment Length'] > 300))]
# add backup check
nifE_backup = nif[(nif.backup_match == 'nifE') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
nifE_backup = nifE_backup.groupby(['GenomeID', 'Hit']).first()
nifE = pd.concat([nifE, nifE_backup])

# only index OG matches
nifN = nif[((nif.residue_match == 'nifN') & (nif.Gene == 'nifN') & (nif.backup_match != 'nifN') & (nif['Alignment Length'] > 300))]
# add backup check
nifN_backup = nif[(nif.backup_match == 'nifN') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
nifN_backup = nifN_backup.groupby(['GenomeID', 'Hit']).first()
nifN = pd.concat([nifN, nifN_backup])

nif = pd.concat([nifH, nifD, nifK, nifE, nifN])
nif.sort_index(inplace = True)

nif.to_csv(f'../results/nif_rescheck_nofilt_{date}.csv')

# Are any nifD,E,K being missed and printed nifN (for example?) should I align all DKEN first and then check for conserved residues?

from cluster_pos import cluster_pos
print('clustering nitrogenase genes', flush=True)

# make sure gene clusters have at least nifHDK
def gene_check(genes):
    if genes.__contains__('nifH'):
        if genes.__contains__('nifD'):
            if genes.__contains__('nifK'):
                return True

# multi-index to cluster by genome, contig
nif = pd.read_csv(f'../results/nif_rescheck_nofilt_{date}.csv')
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

# # filter for genomes to keep
filtered_nif2.set_index(['Hit'], append = True, inplace = True)
filtered_nif2 = filtered_nif2.loc[genomes_to_keep]
filtered_nif2.sort_index(inplace = True)

#clean up cols
filtered_nif2['Gene'] = filtered_nif2['residue_match']
filtered_nif2 = filtered_nif2[['Gene', 'E-value', 'Bit Score', 'Location', 'Orientation', 'Alignment Length', 'Sequence Length', 'GTDB']]
filtered_nif2.drop_duplicates(inplace = True)

# export 
filtered_nif2.to_feather(f'../results/nif_final_{date}.feather')
filtered_nif2.to_csv(f'../results/nif_final_{date}.csv')

nif = pd.read_feather(f'../results/nif_final_{date}.feather')

# export fasta files
nifH = nif[(nif.Gene == 'nifH')]
nifD = nif[(nif.Gene == 'nifD')]
nifK = nif[(nif.Gene == 'nifK')]
nifB = nif[(nif.Gene == 'nifB')]
nifE = nif[(nif.Gene == 'nifE')]
nifN = nif[(nif.Gene == 'nifN')]

gene_list = [nifH, nifD, nifK, nifB, nifE, nifN]
gene_names = ['nifH', 'nifD', 'nifK', 'nifB', 'nifE', 'nifN']

# get fasta sequences for each gene & export to fasta
for gene, name in zip(gene_list, gene_names):
    print(name, flush=True)
    records = []

    for genome,hit in gene.iterrows():

        file = glob.glob(f"../all_rep_proteins_aa/*/{genome[0]}_protein.faa")[0]

        for result in SeqIO.parse(file, "fasta"):
            if result.id == genome[2]:
                # store seq
                gene.loc[genome, 'Seq'] = str(result.seq)
                # convert to seqrecord
                record = SeqRecord(Seq(result.seq), id=genome[0], description=genome[2])
                records.append(record)
                # exit loop once sequence is found
                break
        
    # Write the records to a FASTA file
    with open("../results/check_" + name + ".fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

print("done", flush=True)
