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

nif = pd.read_feather('../results/nif_final_04292025.feather')

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
    print(name)
    records = []

    for genome,hit in gene.iterrows():
        file = glob.glob(f"../all_rep_proteins_aa/*/{genome[0]}_protein.faa")[0]
        for result in SeqIO.parse(file, "fasta"):
            if result.id == genome[-1]:
                # store seq
                gene.loc[genome, 'Seq'] = str(result.seq)
                # convert to seqrecord
                record = SeqRecord(result.seq, id=genome[0], description=genome[-1])
                records.append(record)
                # exit loop once sequence is found
                break
    print(len(records))    
    # Write the records to a FASTA file
    with open("../results/checked_" + name + "04292025.fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
