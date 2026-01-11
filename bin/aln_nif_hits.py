# align hits
from Bio import AlignIO
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
import os
import glob

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np

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

## get fasta sequences for each gene & export to fasta
#for gene, name in zip(gene_list, gene_names):
#    print('getting ' + name + ' fasta', flush=True)
#    records = []
#    for genome,hit in gene.iterrows():
#        hit = hit.Hit
#        file = glob.glob(f"../all_rep_proteins_aa/*/{genome}_protein.faa")[0]
#        
#        for result in SeqIO.parse(file, "fasta"):
#            if result.id == hit:
#                # store seq
#                gene.loc[genome, 'Seq'] = str(result.seq)
#                # convert to seqrecord
#                record = SeqRecord(result.seq, id=genome, description=hit)
#                records.append(record)
#                # exit loop once sequence is found
#                break
        
    # Write the records to a FASTA file
#    with open("../results/" + name + ".fasta", "w") as output_handle:
#        SeqIO.write(records, output_handle, "fasta")

# align fasta files
for gene in gene_names:
    print("aligning "+ gene, flush=True)
    num = eval(f"int({gene}.shape[0]/4000)+1") # how many splits
    os.system(f"seqtk split -n {num} ../results/{gene}_split ../results/{gene}.fasta") # split fasta file
    for i in range(num):
        print(i+1)
        if i+1 < 10:
            os.system(f"seqtk subseq ../results/{gene}.fasta ../results/ref_seq.ids >> ../results/{gene}_split.0000{i+1}.fa") # add reference sequences
            os.system(f"mafft --auto --quiet --thread 4 ../results/{gene}_split.0000{i+1}.fa > ../results/{gene}_split.0000{i+1}.aln")
        else:
            os.system(f"seqtk subseq ../results/{gene}.fasta ../results/ref_seq.ids >> ../results/{gene}_split.000{i+1}.fa") # add reference sequences
            os.system(f"mafft --auto --quiet --thread 4 ../results/{gene}_split.000{i+1}.fa > ../results/{gene}_split.000{i+1}.aln")
