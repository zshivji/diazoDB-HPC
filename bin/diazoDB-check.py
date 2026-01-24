from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import os
import glob
import sys

# reload fasta files?
date = sys.argv[1]

# grab file nif hits
nif_hits = pd.read_csv(f'results/nif_final_{date}.csv')

# read test file
test_file = pd.read_csv('bin/diazoDB-checks.csv')

# compare outputs
for index,row in test_file.iterrows():

    # get test values
    genome = row['genomeID']
    taxonomy = row['taxonomy']
    hit = row['hit']
    gene = row['gene']
    
    # get actual hit from nif_hits
    try:
        actual_hit = nif_hits[(nif_hits['GenomeID'] == genome) & (nif_hits['Hit'] == hit)].Gene.values[0]
        if actual_hit == gene:
            print('PASS: ' + genome + ' ' + taxonomy + ' ' + gene + ' ' + hit, flush=True)
        else:
            print('FAIL: ' + genome + ' ' + taxonomy + ' ' + gene + ' ' + hit, flush=True)
    except:
        print('need to update: ' + genome + ' ' + taxonomy)