import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from Bio import Phylo
import os
import re
import glob
import regex

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from multiprocessing import Pool
from functools import partial
import numpy as np

nif = pd.read_feather('../results/nif_final.feather')

# export fasta files (use .copy() to avoid SettingWithCopyWarning)
nifH = nif[(nif.Gene == 'nifH')].copy()
nifD = nif[(nif.Gene == 'nifD')].copy()
nifK = nif[(nif.Gene == 'nifK')].copy()
nifB = nif[(nif.Gene == 'nifB')].copy()
nifE = nif[(nif.Gene == 'nifE')].copy()
nifN = nif[(nif.Gene == 'nifN')].copy()

gene_list = [nifH, nifD, nifK, nifB, nifE, nifN]
gene_names = ['nifH', 'nifD', 'nifK', 'nifB', 'nifE', 'nifN']

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