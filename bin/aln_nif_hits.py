import pandas as pd
import os
import glob
from tqdm import tqdm
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# get args (reload_fasta)
parser = argparse.ArgumentParser(description='Align Nif gene hits')
parser.add_argument('--reload_fasta', action='store_true', help='Whether to reload fasta sequences (True/False)')
args = parser.parse_args()

# grab both archaea + bacteria hits
nif_archaea = pd.read_feather('../results/archaea/nif.feather')
nif_bacteria = pd.read_feather('../results/bacteria/nif.feather')
nif = pd.concat([nif_archaea, nif_bacteria])

nif.reset_index(inplace = True)
nif.set_index(['GenomeID'], inplace = True)
nif['Seq'] = ''

# separate by annotation
nifH = nif[nif.Gene == 'nifH']
nifD = nif[nif.Gene == 'nifD']
nifK = nif[nif.Gene == 'nifK']
nifB = nif[nif.Gene == 'nifB']
nifE = nif[nif.Gene == 'nifE']
nifN = nif[nif.Gene == 'nifN']

gene_list = [nifH, nifD, nifK, nifB, nifE, nifN]
gene_names = ['nifH', 'nifD', 'nifK', 'nifB', 'nifE', 'nifN']

if args.reload_fasta: # should skip if fasta sequences have already been extracted
    
    # get fasta sequences for each gene & export to fasta
    for gene, name in zip(gene_list, gene_names):
        print('getting ' + name + ' fasta', flush=True)

#        for _ in tqdm(range(len(gene))): # progress bar
        records = []
        for genome,hit in gene.iterrows():
            hit = hit.Hit
            file = glob.glob(f"../all_rep_proteins_aa/*/{genome}_protein.faa")[0]
                
            for result in SeqIO.parse(file, "fasta"):
                if result.id == hit:
                    # store seq
                    gene.loc[genome, 'Seq'] = str(result.seq)
                    # convert to seqrecord
                    record = SeqRecord(result.seq, id=genome, description=hit)
                    records.append(record)
                    # exit loop once sequence is found
                    break
                    
            # Write the records to a FASTA file
        with open("../results/intermediate_fastas/" + name + ".fasta", "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

#align fasta files
for gene in gene_names:
    print("aligning "+ gene, flush=True)
    num = eval(f"int({gene}.shape[0]/200)+1") # how many splits
    os.system(f"seqtk split -n {num} ../results/fasta_splits/{gene}_split ../results/intermediate_fastas/{gene}.fasta") # split fasta file
    for i in range(num):
        os.system(f"seqtk subseq ../results/intermediate_fastas/{gene}.fasta ../results/ref_seq.ids >> ../results/fasta_splits/{gene}_split.{str(i+1).zfill(5)}.fa") # add reference sequences
        os.system(f"mafft --auto --quiet --thread 4 ../results/fasta_splits/{gene}_split.{str(i+1).zfill(5)}.fa > ../results/fasta_splits/{gene}_split.{str(i+1).zfill(5)}.aln")
