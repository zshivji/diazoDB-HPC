import pandas as pd
import os
import glob
from tqdm import tqdm
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from helper import load_nif

# get args (reload_fasta)
parser = argparse.ArgumentParser(description='Align Nif gene hits')
parser.add_argument('--reload_fasta', action='store_true', help='Whether to reload fasta sequences (True/False)')
args = parser.parse_args()

# grab both archaea + bacteria hits
hits = load_nif(update_index = ['GenomeID']) # load nif dataframe
hits['Seq'] = ''

# separate by annotation
nifH = hits[hits.Gene == 'nifH']
nifD = hits[hits.Gene == 'nifD']
nifK = hits[hits.Gene == 'nifK']
nifB = hits[hits.Gene == 'nifB']
nifE = hits[hits.Gene == 'nifE']
nifN = hits[hits.Gene == 'nifN']
anfO = hits[hits.Gene == 'anfO']
nifG = hits[hits.Gene == 'anfG' or hits.Gene == 'vnfG']

gene_list = [nifH, nifD, nifK, nifB, nifE, nifN, anfO, nifG]
gene_names = ['nifH', 'nifD', 'nifK', 'nifB', 'nifE', 'nifN', 'anfO', 'nifG']

if args.reload_fasta: # should skip if fasta sequences have already been extracted
    
    # get fasta sequences for each gene & export to fasta
    for gene, name in zip(gene_list, gene_names):
        print('getting ' + name + ' fasta', flush=True)

#        for _ in tqdm(range(len(gene))): # progress bar
        records = []
        for genome,hit in gene.iterrows():
            hit = hit.Hit
            file = glob.glob(f"../protein_faa_reps/*/{genome}_protein.faa")[0]
                
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
        os.makedirs("../results/tmp", exist_ok=True)
        with open("../results/tmp/" + name + ".fasta", "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

#align fasta files
for gene in gene_names:
    print("aligning "+ gene, flush=True)
    num = eval(f"int({gene}.shape[0]/200)+1") # how many splits
    os.makedirs("../results/tmp/fasta_splits", exist_ok=True)
    os.system(f"seqtk split -n {num} ../results/tmp/fasta_splits/{gene}_split ../results/tmp/{gene}.fasta") # split fasta file
    for i in range(num):
        os.system(f"seqtk subseq ../results/tmp/{gene}.fasta ../results/ref_seq.ids >> ../results/tmp/fasta_splits/{gene}_split.{str(i+1).zfill(5)}.fa") # add reference sequences
        os.system(f"mafft --auto --quiet --thread 4 ../results/tmp/fasta_splits/{gene}_split.{str(i+1).zfill(5)}.fa > ../results/tmp/fasta_splits/{gene}_split.{str(i+1).zfill(5)}.aln")
