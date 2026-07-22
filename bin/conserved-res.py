import pandas as pd
import numpy as np
import glob
import json
import os
from tqdm import tqdm
import argparse

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from helper import load_nif

# get args (reload_fasta)
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make and align temporary fasta files based on filtered hits."
            "Check for conserved residues in alignments."
            "Update nif gene assignments based on conserved residue matching."
        )
    )
    parser.add_argument(
        "--reload_fasta",
        help="Whether or not to reload the fasta files (slow but necessary if hmmsearch/parse_hmmp.py changed).",
        action='store_true'
    )

    parser.add_argument(
        "--align_fasta",
        help="Whether or not to align the fasta files (slow but necessary if hmmsearch/parse_hmmp.py changed). Automatically defaults to True if --reload_fasta is specified.",
        action='store_true'
    )

    parser.add_argument(
        "--config_file",
        help="Path to the JSON configuration file containing reference sequences and important residues. Default bin/nif-config.json",
        default='nif-config.json'
    )

    return parser.parse_args()

# grab both archaea + bacteria hits
hits = load_nif(update_index = ['GenomeID']) # load nif dataframe
hits['Seq'] = ''

# read json
config = json.load(open('nif-config.json', 'r'))

# separate hits by gene
gene_list = {gene: hits[hits.Gene == gene] for gene in config.keys()}
print(f"checking conserved residues for: {list(gene_list)}", flush=True)

def write_fasta(): # get fasta sequences for each gene & export to fasta        
    for name, gene in gene_list.items():
        print('getting ' + name + ' fasta', flush=True)
        # filter out hits that don't meet length criteria (domain is too short, pre-filtering decreases runtime)
        minl = config[name]['min-length']
        gene = gene[abs(gene['Domain Start'] - gene['Domain End']) >= minl]
        gene_list[name] = gene #### DOESN"T UPDATE GLOBAL DICT, also need to export w/ sequence col so that it runs even if write_fasta skipped

        records = []
        for genome,hit in gene.iterrows():
            file = glob.glob(f"../protein_faa_reps_232/*/{genome}_protein.faa")[0]
                
            for result in SeqIO.parse(file, "fasta"):
                if result.id == hit.Hit[:-2]: # find matching hit (dropping domain suffix)
                    # store seq
                    domain_seq = result.seq[hit['Domain Start']:hit['Domain End']] # trim seq to match domain len
                    gene.loc[genome, 'Seq'] = str(domain_seq)
                    # convert to seqrecord
                    record = SeqRecord(domain_seq, id=genome, description=hit.Hit)
                    records.append(record)
                    # exit loop once sequence is found
                    break
        print(len(records), "records found for " + name, flush=True)            
        # Write the records to a FASTA file
        os.makedirs("../results/tmp", exist_ok=True)
        with open("../results/tmp/" + name + ".fasta", "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

def aln_fasta():
    #align fasta files
    for name, gene in gene_list.items():
        print("aligning "+ name, flush=True)
        num = int(gene.shape[0]/200)+1 # how many splits
        os.makedirs("../results/tmp/fasta_splits", exist_ok=True)
        os.system(f"seqtk split -n {num} ../results/tmp/fasta_splits/{name}_split ../results/tmp/{name}.fasta") # split fasta file
        for i in range(num):
            os.system(f"seqtk subseq ../results/tmp/{name}.fasta ../results/ref_seq.ids >> ../results/tmp/fasta_splits/{name}_split.{str(i+1).zfill(5)}.fa") # add reference sequences
            os.system(f"mafft --auto --quiet --thread 4 ../results/tmp/fasta_splits/{name}_split.{str(i+1).zfill(5)}.fa > ../results/tmp/fasta_splits/{name}_split.{str(i+1).zfill(5)}.aln")


def check_gene(file, ref_seq, important_residues, residue_scores, passing_score, p=False):
    if len(important_residues) != len(residue_scores):
        raise ValueError(
            f"{file[-20:-16]}: residues and residue_scores must have the same length "
            f"({len(important_residues)} != {len(residue_scores)})"
        )

    residue_to_score = dict(zip(important_residues, residue_scores))
    alignment = AlignIO.read(f"{file}", "fasta")

    for record in alignment:
        if ref_seq in record.description:
            aln = record.seq
            break

    # map important cols in aln to ref_seq
    residue_to_alignment = {}
    residue_idx = 0  # index in original (ungapped) sequence
    col2res = {}
    col2score = {}

    for aln_idx, char in enumerate(aln):
        if char != '-':
            residue_idx += 1
            if residue_idx in important_residues:
                residue_to_alignment[residue_idx] = aln_idx
                col2res[aln_idx] = aln[aln_idx]
                col2score[aln_idx] = residue_to_score[residue_idx]
                
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
                score += 1 + col2score[col]
            elif row[col] != '-': # higher weight for correct C
                    score += 1
        if score >= passing_score:
            return score
        else:
            return np.nan

    pssm['score'] = pssm.apply(check_res, axis = 1)
    pssm.dropna(subset = ['score'], inplace = True)
    pssm.rename(columns=col2res, inplace=True)

    return pssm

def check_gene_wrapper():
    print("Checking for conserved residues...", flush=True)

    # performed conserved residue matching
    all_genes_checked = {}

    for name, gene in gene_list.items():
        print(f'checking {name}', flush=True)
        
        # clear df
        checked = pd.DataFrame()

        # load json parameters
        ref_seq = config[name]['ref_gene']
        important_residues = config[name]['residues']
        residue_scores = config[name].get('residue_scores', config[name].get('reside_scores'))
        if residue_scores is None:
            residue_scores = [1] * len(important_residues)
        passing_score = config[name]['passing_score']

        # initialize dataframe
        checked = check_gene(f'../results/tmp/fasta_splits/{name}_split.00001.aln', ref_seq, important_residues, residue_scores, passing_score, p=True)

        for file in glob.glob(f'../results/tmp/fasta_splits/{name}_split.*.aln'):
            if f'{name}_split.00001' in file:
                continue
            new = check_gene(file, ref_seq, important_residues, residue_scores, passing_score)
            checked = pd.concat([checked, new])

        checked.drop_duplicates(subset = ['hit'], inplace = True)
        checked.set_index(['hit'], append= True, inplace = True)

        # potentially add code for (1) keeping best hit per genome and (2) removing hits that are in other genes (e.g. remove nifD,E hits from nifK)

        all_genes_checked[name] = checked

        print(str(len(checked)) + f" {name} seqs")

    return all_genes_checked

# #### EXPORT

def export_results(all_genes_checked):

    # reload nif dataframe
    nif = load_nif(update_index = ['GenomeID', 'Hit','Gene']) # load nif dataframe

    # append updated annotation (based on conserved residue matching) to nif files
    nif['residue_match'] = ''

    nif.reset_index(level='Gene', inplace = True)
    nif.sort_index(inplace = True)

    # update residue match column in nif df
    for name, gene in gene_list.items():
        for genome, cols in eval(f"all_genes_checked[nif{gene}].iterrows()"):
            nif.loc[(genome[0], genome[1]), 'residue_match'] = "nif" + gene
        
        # filter to get hits that passed residue matching
        gene = f"nif{gene}"
        maxl = config[gene]['max-length']
        minl = config[gene]['min-length']
        exec(f"{gene} = nif[((nif.residue_match == '{gene}') & (nif.Gene == '{gene}') & (nif['Alignment Length'] >= {minl}) & (nif['Alignment Length'] < {maxl}))]", globals())

    nif = pd.concat([nifH, nifD, nifK, nifE, nifN, nifNB, nifB, anfO, nifG])
    nif.sort_index(inplace = True)

    nif.to_csv(f'../results/tmp/nif_rescheck_nofilt.csv')

#### END OF CONSERVED RESIDUE CHECKING, BELOW IS VALIDATion/BACKUP CODE ####





# #backup check --> likely don't need this since we're reducing alignments to 200 seq per file
# #get seq that failed checks
# print('getting failed sequences', flush=True)

# seqs = []

# for gene in 'DKEN':

#     # for each gene, get all hmm hit acc
#     result = list(SeqIO.parse(f"../results/fasta_splits/nif{gene}.fasta", "fasta"))
#     hit = [record.description.split(" ")[-1] for record in result]

#     # get seq that failed check
#     for record, acc in zip(result, hit):
#         if acc not in eval(f"list(nif{gene}_checked.index.get_level_values(1).unique())"):
#             seq = SeqRecord(Seq(record.seq), id=record.id, description=acc)
#             seqs.append(record)

# print(str(len(seqs)) + " seqs failed nifDKEN checks", flush=True) # counting how many failed pre/post 2026 edits


def main() -> None:
    args = parse_args()

    # reload + align fasta files if specified
    if args.reload_fasta:
        write_fasta()
        args.align_fasta = True  # Automatically align if reloading fasta
    if args.align_fasta:
        aln_fasta()
    # check for conserved residues
    checked = check_gene_wrapper()
    # export
    export_results(checked)

if __name__ == "__main__":
    main()

# # Write the records to a FASTA file
# with open(f"../results/fasta_splits/nifDKEN.fasta", "w") as output_handle:
#     SeqIO.write(seqs, output_handle, "fasta")

# # add reference sequences
# for gene in 'DKEN':
#     os.system(f"seqtk subseq ../results/fasta_splits/nif{gene}.fasta ../results/ref_seq.ids >> ../results/fasta_splits/nifDKEN.fasta")

# print('aligning failed sequences', flush=True)
# # aln all seqs
# num = int(len(seqs)/200) +1  # how many splits
# os.system(f"seqtk split -n {num} ../results/fasta_splits/nifDKEN_split ../results/fasta_splits/nifDKEN.fasta") # split fasta file
# for i in range(num):
#     print(i+1)
#     os.system(f"seqtk subseq ../results/fasta_splits/nifDKEN.fasta ../results/ref_seq.ids >> ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.fa") # add reference sequences
#     os.system(f"mafft --auto --quiet --thread 4 ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.fa > ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.aln")
    

# backup check --> likely don't need this since we're reducing alignments to 200 seq per file
# print('checking failed sequences', flush=True)

# nifD_backup = check_gene('nifDKEN_split.00001', ref_seq_nifD, important_residues_nifD, passing_score)
# nifK_backup = check_gene('nifDKEN_split.00001', ref_seq_nifK, important_residues_nifK, passing_score)
# nifE_backup = check_gene('nifDKEN_split.00001', ref_seq_nifE, important_residues_nifE, passing_score)
# nifN_backup = check_gene('nifDKEN_split.00001', ref_seq_nifN, important_residues_nifN, passing_score)

# for file in glob.glob(f'nifDKEN_split.00*.aln'):
#     if file == '../results/fasta_splits/nifDKEN_split.00001.aln':
#         continue
#     nifD_backup = check_gene(file[:-4], ref_seq_nifD, important_residues_nifD, passing_score)
#     nifD_backup = pd.concat([nifD_backup, new])

#     nifK_backup = check_gene(file[:-4], ref_seq_nifK, important_residues_nifK, passing_score)
#     nifK_backup = pd.concat([nifK_backup, new])

#     nifE_backup = check_gene(file[:-4], ref_seq_nifE, important_residues_nifE, passing_score)
#     nifE_backup = pd.concat([nifE_backup, new])

#     nifN_backup = check_gene(file[:-4], ref_seq_nifN, important_residues_nifN, passing_score)
#     nifN_backup = pd.concat([nifN_backup, new])

# # set index as genome, hit
# nifD_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per contig
# nifD_backup['gene_cluster'] = 0
# nifD_backup = nifD_backup.loc[nifD_backup.groupby(['contig'])['score'].idxmax()]
# nifD_backup.drop_duplicates(inplace = True)

# print(str(len(nifD_backup.index.unique())) + " nifD seqs", flush=True)
# nifD_checked = pd.concat([nifD_checked, nifD_backup])

# # set index as genome, hit
# nifK_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per contig
# nifK_backup['gene_cluster'] = 0
# nifK_backup = nifK_backup.loc[nifK_backup.groupby(['contig'])['score'].idxmax()]
# nifK_backup.drop_duplicates(inplace = True)

# print(str(len(nifK_backup.index.unique())) + " nifK seqs", flush=True)
# nifK_checked = pd.concat([nifK_checked, nifK_backup])

# # set index as genome, hit
# nifE_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per gene cluster
# nifE_backup['gene_cluster'] = 0
# nifE_backup = nifE_backup.loc[nifE_backup.groupby(['contig'])['score'].idxmax()]
# nifE_backup.drop_duplicates(inplace = True)

# print(str(len(nifE_backup.index.unique())) + " nifE seqs", flush=True)
# nifE_checked = pd.concat([nifE_checked, nifE_backup])

# # set index as genome, hit
# nifN_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per gene cluster
# nifN_backup['gene_cluster'] = 0
# nifN_backup = nifN_backup.loc[nifN_backup.groupby(['contig'])['score'].idxmax()]
# nifN_backup.drop_duplicates(inplace = True)

# print(str(len(nifN_backup.index.unique())) + " nifN seqs", flush=True)
# nifN_checked = pd.concat([nifN_checked, nifN_backup])

# # append updated annotation (based on conserved residue matching) to nif files
# nif['residue_match'] = ''
# nif['backup_match'] = ''

# nif.set_index(['Hit'], append = True, inplace = True)
# nif.sort_index(inplace = True)

# # update residue match column in nif df
# for gene in 'HDKEN':
#     for genome, cols in eval(f"nif{gene}_checked.iterrows()"):
#         nif.loc[(genome, cols.hit), 'residue_match'] = "nif" + gene

#  # add backup check       
# for gene in 'DKEN':
#     for genome, cols in eval(f"nif{gene}_backup.iterrows()"):
#         nif.loc[(genome[0], genome[1]), 'backup_match'] = "nif" + gene

# # filter to get hits that passed residue matching
# nifH = nif[(nif.residue_match == 'nifH') & (nif.Gene == 'nifH') & (nif['Alignment Length'] > 200)]
# # nifB = nif[(nif.Gene == 'nifB')] # not done
# # nifB['residue_match'] = 'nifB'

# # only index OG matches
# nifD = nif[((nif.residue_match == 'nifD') & (nif.Gene == 'nifD') & (nif.backup_match != 'nifD') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifD_backup = nif[(nif.backup_match == 'nifD') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifD_backup = nifD_backup.groupby(['GenomeID', 'Hit']).first()
# nifD = pd.concat([nifD, nifD_backup])

# # only index OG matches
# nifK = nif[((nif.residue_match == 'nifK') & (nif.Gene == 'nifK') & (nif.backup_match != 'nifK') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifK_backup = nif[(nif.backup_match == 'nifK') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifK_backup = nifK_backup.groupby(['GenomeID', 'Hit']).first()
# nifK = pd.concat([nifK, nifK_backup])

# # only index OG matches
# nifE = nif[((nif.residue_match == 'nifE') & (nif.Gene == 'nifE') & (nif.backup_match != 'nifE') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifE_backup = nif[(nif.backup_match == 'nifE') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifE_backup = nifE_backup.groupby(['GenomeID', 'Hit']).first()
# nifE = pd.concat([nifE, nifE_backup])

# # only index OG matches
# nifN = nif[((nif.residue_match == 'nifN') & (nif.Gene == 'nifN') & (nif.backup_match != 'nifN') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifN_backup = nif[(nif.backup_match == 'nifN') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifN_backup = nifN_backup.groupby(['GenomeID', 'Hit']).first()
# nifN = pd.concat([nifN, nifN_backup])

# nif = pd.concat([nifH, nifD, nifK, nifE, nifN])
# nif.sort_index(inplace = True)

# nif.to_csv(f'../results/fasta_splits/nif_rescheck_nofilt.csv')

# # Are any nifD,E,K being missed and printed nifN (for example?) should I align all DKEN first and then check for conserved residues?

# EXPORT
