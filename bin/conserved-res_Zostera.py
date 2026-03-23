from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import numpy as np
import os
import glob
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

dir = sys.argv[1]
date = sys.argv[2]

def check_gene(gene, ref_seq, important_residues, passing_score, p=False):

    alignment = AlignIO.read(f"{gene}.aln", "fasta")

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
            print(f"Residue {residue} corresponds to alignment index {aln_idx}: {aln[aln_idx]}")

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
    # rename cols
    pssm.rename(columns = col2res, inplace = True)

    return pssm

print("running", flush=True)

# arsA
print("checking arsA", flush=True)

important_residues_arsA = [113, 172, 422] # CCC in e. coli

ref_seq_arsA = 'tr|A0A2A5MC03|A0A2A5MC03_9ENTR'

passing_score = 9

# initialize dataframe
arsA_checked = check_gene(f'../results/fasta_splits/arsA_{dir}_split.00001', ref_seq_arsA, important_residues_arsA, passing_score, p=True)

for file in glob.glob(f'../results/fasta_splits/arsA_{dir}_split.00*.aln'):
    if file == f'../results/fasta_splits/arsA_{dir}_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_arsA, important_residues_arsA, passing_score)
    arsA_checked = pd.concat([arsA_checked, new])

arsA_checked.drop_duplicates(subset = ['hit'], inplace = True)
arsA_checked.set_index(['hit'], append= True, inplace = True)

print(str(len(arsA_checked)) + " arsA seqs")

# append updated annotation (based on conserved residue matching) to ars files

# grab both archaea + bacteria hits
if dir == 'WW':
    ars = pd.read_feather(f'../results/{dir}/tophits.feather')
else:
    ars_archaea = pd.read_feather(f'../results/archaea/tophits.feather')
    ars_bacteria = pd.read_feather('../results/bacteria/tophits.feather')
    ars = pd.concat([ars_archaea, ars_bacteria])

ars.reset_index(inplace = True)
ars.set_index(['GenomeID', 'Hit'], inplace = True)
ars['Seq'] = ''
ars['residue_match'] = ''
ars['backup_match'] = ''

ars.sort_index(inplace = True)

# update residue match column in ars df
for genome, cols in eval(f"arsA_checked.iterrows()"):
    ars.loc[(genome[0], genome[1]), 'residue_match'] = "arsA"

# filter to get hits that passed residue matching
ars.reset_index(level = 'Hit', inplace = True)
ars.set_index('contig', append=True, inplace = True)

arsA = ars[(ars.residue_match == 'arsA') & (ars.Gene == 'arsA') & (ars['Alignment Length'] > 400)]
ars = ars[ars.Gene != 'arsA'] # remove all arsA hits

arsB = ars[(ars.Gene == 'arsB') & (ars['Alignment Length'] >= 400) & (ars['Alignment Length'] <= 500)]
ars = ars[ars.Gene != 'arsB'] # remove all arsB hits

ars = pd.concat([arsA, arsB, ars]) # add checked arsA hits back
ars.sort_index(inplace = True)

ars.to_csv(f'../results/ars_{dir}_rescheck_nofilt_{date}.csv', index = True)

from cluster_pos import cluster_pos

# make sure gene clusters have at least arsRBC
def gene_check(genes):
    if genes.__contains__('arsR'):
        if genes.__contains__('arsB'):
            if genes.__contains__('arsC'):
                return True

## multi-index to cluster by genome, contig
ars = pd.read_csv(f'../results/ars_{dir}_rescheck_nofilt_{date}.csv', index_col=['GenomeID', 'contig'])
# drop duplicate hits, but keep most sig e-value
ars.sort_values(by=['E-value'], inplace=True)
ars.drop_duplicates(inplace = True, subset=['Hit'], keep='first')
ars.sort_index(inplace = True) # improve performance

# filter for genome, contig with at least 3 unique genes (arsRBC)
filtered_ars = ars.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Gene'].nunique() >= 3)

genomes_to_keep = []
# iterate through each genome and contig
for genome in filtered_ars.index.get_level_values(0).unique(): # iterate through each genome
    for contig in filtered_ars.loc[genome].index.get_level_values(0).unique(): # iterate through each contig

        tmp = filtered_ars.loc[(genome, contig)]

        # only keep numbers that have clusters >= 3
        pos_clusters = cluster_pos(tmp.Hit.unique(), 10)

        # save clusters that have at least arsRBC
        for ind, cl in enumerate(pos_clusters):
            pos = [contig + '_' + str(p) for p in cl]
            no_pos = len(pos)
            
            # need at least 3 genes to continue
            if no_pos < 3:
                continue

            # only keep hits that are in the cluster
            tmp2 = tmp[tmp.Hit.isin(pos)].reset_index()

            # check if all genes are present
            if gene_check(tmp2.Gene.to_list()):
                # get index
                items = [(genome, contig, hit) for hit in tmp2.Hit]
                genomes_to_keep.extend(items)

# filter for genomes to keep
filtered_ars.set_index(['Hit'], append = True, inplace = True)
filtered_ars = filtered_ars.loc[genomes_to_keep]
filtered_ars.sort_index(inplace = True)

#clean up cols
filtered_ars = filtered_ars[['Gene', 'E-value', 'Bit Score', 'Location', 'Orientation', 'Alignment Length', 'Sequence Length', 'GTDB']]
filtered_ars.drop_duplicates(inplace = True)

# export 
filtered_ars.to_feather(f'../results/ars_{dir}_final_{date}.feather')
filtered_ars.to_csv(f'../results/ars_{dir}_final_{date}.csv')
