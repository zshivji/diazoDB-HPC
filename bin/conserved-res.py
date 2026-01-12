from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import numpy as np
import os
import glob

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

    # check if cols contain correct residue for function
    def check_res(row):
        score = 0
        for col in residue_to_alignment.values():
            if row[col] == aln[col]:
                score += 2
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

print("running", flush=True)
nif_archaea = pd.read_feather('../results/archaea/nif.feather')
nif_bacteria = pd.read_feather('../results/bacteria/nif.feather')
nif = pd.concat([nif_archaea, nif_bacteria])

nif.reset_index(inplace = True)
nif.set_index(['GenomeID'], inplace = True)
nif['Seq'] = ''

# nifH
print("checking nifH", flush=True)
#important_residues_nifH = [8, 9, 10, 11, 13, 14, 15, 38, 40, 42, 95, 127, 130] # methanosarcina
important_residues_nifH = [10, 11, 12, 13, 15, 16, 17, 40, 42, 44, 98, 130, 133] # azotobacter

#ref_seq = 'NC_003552.1_4100'
ref_seq_nifH = 'NZ_BSFG01000004.1_46'

passing_score = 24

nifH_checked = check_gene('nifH_split.00001', ref_seq_nifH, important_residues_nifH, passing_score, p=True)

for file in glob.glob(f'../results/nifH_split.000*.aln'):
    if file == '../results/nifH_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifH, important_residues_nifH, passing_score)
    nifH_checked = pd.concat([nifH_checked, new])

print(str(len(nifH_checked)) + " nifH seqs")

# nifD
print('checking nifD', flush=True)
# grab aln numbering to map to ref seq (azoto, methano)
#important_residues_nifD = [58, 84, 149, 186, 190, 267, 485] # methanosarcina
important_residues_nifD = [62, 88, 154, 191, 195, 275, 442] # azotobacter

#ref_seq = 'NC_003552.1_4103'
ref_seq_nifD = 'NZ_BSFG01000004.1_45' # azo

passing_score = 13 # appears that some groups have H442R,K

nifD_checked = check_gene('nifD_split.00001', ref_seq_nifD, important_residues_nifD, passing_score, p=True)

for file in glob.glob(f'../results/nifD_split.000*.aln'):
    if file == '../results/nifD_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifD, important_residues_nifD, passing_score)
    nifD_checked = pd.concat([nifD_checked, new])

print(str(len(nifD_checked)) + " nifD seqs", flush=True)

# nifE 
print('checking nifE', flush=True)
# grab aln numbering to map to ref seq (azoto, methano)
# important_residues_nifE = [34, 46, 71, 132] # methanosarcina
important_residues_nifE = [25, 37, 62, 124] # azo, known catalytic
#important_residues_nifE = [25, 37, 62, 124, 46, 51, 155, 163, 324, 395]
# important_residues = [25, 37, 62, 124] # azotobacter YF39, G31, G41, L46, D51, E90, G96, L102, P115, I129, D132, P147, G155, G163, G202, G221, P268, G276, G324, G347, D395, G427

# ref_seq = 'NC_003552.1_4105' # methano
ref_seq_nifE = 'NZ_BSFG01000004.1_39' # azo

passing_score = 8 # could be more rigorous

nifE_checked = check_gene('nifE_split.00001', ref_seq_nifE, important_residues_nifE, passing_score, p=True)

for file in glob.glob(f'../results/nifE_split.000*.aln'):
    if file == '../results/nifE_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifE, important_residues_nifE, passing_score)
    nifE_checked = pd.concat([nifE_checked, new])

print(str(len(nifE_checked)) + " nifE seqs", flush=True)

# nifK
print('checking nifK', flush=True)
# grab aln numbering to map to ref seq (azoto, methano)
#important_residues = [23, 48, 106, 141] # methanosarcina
important_residues_nifK = [70, 95, 153, 188] # azotobacter

#ref_seq = 'NC_003552.1_4104'
ref_seq_nifK = 'NZ_BSFG01000004.1_44'

passing_score = 7 # S188 has not been verified to be conserved (8 is too conservative)

nifK_checked = check_gene('nifK_split.00001', ref_seq_nifK, important_residues_nifK, passing_score, p=True)

for file in glob.glob(f'../results/nifK_split.000*.aln'):
    if file == '../results/nifK_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifK, important_residues_nifK, passing_score)
    nifK_checked = pd.concat([nifK_checked, new])

print(str(len(nifK_checked)) + " nifK seqs", flush=True)

# nifN --> any other residues??
print('checking nifN', flush=True)
# grab aln numbering to map to ref seq (azoto, methano)
#important_residues = [] # methanosarcina
important_residues_nifN = [44] # azotobacter

#ref_seq = 'NC_003552.1_4106'
ref_seq_nifN = 'NZ_BSFG01000004.1_38'

passing_score = 2

nifN_checked = check_gene('nifN_split.00001', ref_seq_nifN, important_residues_nifN, passing_score, p=True)

for file in glob.glob(f'../results/nifN_split.000*.aln'):
    if file == '../results/nifN_split.00001.aln':
        continue
    new = check_gene(file[:-4], ref_seq_nifN, important_residues_nifN, passing_score)
    nifN_checked = pd.concat([nifN_checked, new])

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
        if acc not in eval(f"list(nif{gene}_checked.hit.unique())"):
            seq = SeqRecord(Seq(record.seq), id=record.id, description=acc)
            seqs.append(record)

print(str(len(seqs)) + " seqs failed nifDKEN check", flush=True)

# Write the records to a FASTA file
with open("../results/nifDKEN.fasta", "w") as output_handle:
    SeqIO.write(seqs, output_handle, "fasta")

# add reference sequences
for gene in 'DKEN':
    os.system(f"seqtk subseq ../results/nif{gene}.fasta ../results/ref_seq.ids >> ../results/nifDKEN.fasta")

print('aligning failed sequences', flush=True)
# aln all seqs
num = int(len(seqs)/4000) +1  # how many splits
os.system(f"seqtk split -n {num} ../results/nifDKEN_split ../results/nifDKEN.fasta") # split fasta file
for i in range(num):
    if i+1 < 10:
        os.system(f"seqtk subseq ../results/nifDKEN.fasta ../results/ref_seq.ids >> ../results/nifDKEN_split.0000{i+1}.fa") # add reference sequences
        os.system(f"mafft --auto --quiet --thread 4 ../results/nifDKEN_split.0000{i+1}.fa > ../results/nifDKEN_split.0000{i+1}.aln")
    else:
        os.system(f"seqtk subseq ../results/nifDKEN.fasta ../results/ref_seq.ids >> ../results/nifDKEN_split.000{i+1}.fa") # add reference sequences
        os.system(f"mafft --auto --quiet --thread 4 ../results/nifDKEN_split.000{i+1}.fa > ../results/nifDKEN_split.000{i+1}.aln")


# backup check
print('checking failed sequences', flush=True)

nifD_backup = check_gene('nifDKEN_split.00001', ref_seq_nifD, important_residues_nifD, passing_score, p=True)
nifK_backup = check_gene('nifDKEN_split.00001', ref_seq_nifK, important_residues_nifK, passing_score, p=True)
nifE_backup = check_gene('nifDKEN_split.00001', ref_seq_nifE, important_residues_nifE, passing_score, p=True)
nifN_backup = check_gene('nifDKEN_split.00001', ref_seq_nifN, important_residues_nifN, passing_score, p=True)

for file in glob.glob(f'nifDKEN_split.000*.aln'):
    if file == '../results/nifDKEN.00001.aln':
        continue
    nifD_backup = check_gene(file[:-4], ref_seq_nifD, important_residues_nifD, 14)
    nifD_backup = pd.concat([nifD_backup, new])

    nifK_backup = check_gene(file[:-4], ref_seq_nifK, important_residues_nifK, 8)
    nifK_backup = pd.concat([nifK_backup, new])

    nifE_backup = check_gene(file[:-4], ref_seq_nifE, important_residues_nifE, 7)
    nifE_backup = pd.concat([nifE_backup, new])

    nifN_backup = check_gene(file[:-4], ref_seq_nifN, important_residues_nifN, 2)
    nifN_backup = pd.concat([nifN_backup, new])


print(str(len(nifD_backup)) + " nifD seqs", flush=True)
nifD_checked = pd.concat([nifD_checked, nifD_backup])

nifK_backup = nifK_backup[~nifK_backup.hit.isin(nifD_checked.hit.to_list())] # remove nifD
print(str(len(nifK_backup)) + " nifK seqs", flush=True)
nifK_checked = pd.concat([nifK_checked, nifK_backup])

nifE_backup = nifE_backup[~nifE_backup.hit.isin(nifD_checked.hit.to_list())] # remove nifD
nifE_backup = nifE_backup[~nifE_backup.hit.isin(nifK_checked.hit.to_list())] # remove nifK
print(str(len(nifE_backup)) + " nifE seqs", flush=True)
nifE_checked = pd.concat([nifE_checked, nifE_backup])

nifN_backup = nifN_backup[~nifN_backup.hit.isin(nifD_checked.hit.to_list())] # remove nifD
nifN_backup = nifN_backup[~nifN_backup.hit.isin(nifK_checked.hit.to_list())] # remove nifK
nifN_backup = nifN_backup[~nifN_backup.hit.isin(nifE_checked.hit.to_list())] # remove nifK
print(str(len(nifN_backup)) + " nifN seqs", flush=True)
nifN_checked = pd.concat([nifN_checked, nifN_backup])

# append updated annotation (based on conserved residue matching) to nif files
nif['residue_match'] = ''

nif.set_index(['Hit'], append = True, inplace = True)
nif.sort_index(inplace = True)

# update residue match column in nif df
for gene in 'HDKEN':
    for genome, cols in eval(f"nif{gene}_checked.iterrows()"):
        nif.loc[(genome, cols.hit), 'residue_match'] = "nif" + gene

nif.to_csv('../results/nif_residue_match_no_filter.csv')

# make sure original annotation matches residue matching annotation
nifH = nif[(nif.residue_match == 'nifH')]
nifD = nif[(nif.residue_match == 'nifD')]
nifK = nif[(nif.residue_match == 'nifK')]
nifB = nif[(nif.residue_match == 'nifB')] # not done
nifE = nif[(nif.residue_match == 'nifE')]
nifN = nif[(nif.Gene == 'nifN') & (nif.Gene == nif.residue_match)]

nif = pd.concat([nifH, nifD, nifK, nifB, nifE, nifN])

# Are any nifD,E,K being missed and printed nifN (for example?) should I align all DKEN first and then check for conserved residues?

from cluster_pos import cluster_pos
print('clustering nitrogenase genes', flush=True)
# make sure gene clusters have at least nifHDK/nifHEK/nifHDN/nifHEN or more
def gene_check(genes):
    if genes.__contains__('nifH'):
        if genes.__contains__('nifD') or genes.__contains__('nifE'):
            if genes.__contains__('nifK') or genes.__contains__('nifN'):
                return True

# multi-index to cluster by genome, contig
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
        pos_clusters = cluster_pos(tmp.Hit.unique())
        
        # for each cluster, make sure nifHDKEN (nifHDK)
        for cl in pos_clusters:
            pos = [contig + '_' + str(p) for p in cl]
            no_pos = len(pos)

            # only keep hits that are in the cluster
            tmp2 = tmp[tmp.Hit.isin(pos)]

            # check if all genes are present
            if gene_check(tmp2.Gene.to_list()):

                # get index
                items = [(genome, contig, hit) for hit in tmp2.Hit]

                genomes_to_keep.extend(items)

# # filter for genomes to keep
filtered_nif2.set_index(['Hit'], append = True, inplace = True)
filtered_nif2 = filtered_nif2.loc[genomes_to_keep]
filtered_nif2.sort_index(inplace = True)

# export 
filtered_nif2.to_feather('../results/check_nif_04182025.feather')
filtered_nif2.to_csv('../results/check_nif_04182025.csv')

nif = pd.read_feather('../results/check_nif_04182025.feather')

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
