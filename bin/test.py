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

def test(checked: pd.DataFrame) -> None:
    """
    Compare the nif hits in the final nif file with the expected hits in the test file.
    """
    # read test file
    test_file = pd.read_csv('diazoDB-checks.csv')
    test_file.set_index(['genomeID', 'taxonomy'], inplace=True)

    checked.reset_index(inplace=True)

    # create log file
    with open('diazoDB-check-log.txt', 'w') as log:
        # iterate through genomes in test file
        for genome in test_file.index.get_level_values(0).unique():
            taxonomy = test_file.loc[genome].index.get_level_values(0)[0]
            # write genome, taxonomy to log
            log.write(f'Checking {genome} {taxonomy}\n')

            # compare if correct # genes found
            checked_count = checked[checked['GenomeID'] == genome].shape[0]
            test_file_count = test_file.loc[genome].shape[0]

            if checked_count == test_file_count:
                log.write('PASS: correct number of genes found\n')
            elif checked_count < test_file_count:
                log.write(f'FAIL: missing {test_file_count - checked_count} hits \n')
            elif checked_count > test_file_count:
                    extra = [hit for hit in checked[checked['GenomeID'] == genome].protein.to_list() if hit not in test_file.loc[genome].hit.to_list()]
                    log.write(f'FAIL: {checked_count - test_file_count} extra hits: {extra} \n')

            # compare gene annotations
            for index,row in test_file.loc[genome].iterrows():

                # get test values
                hit = row['hit']
                actual_hit = row['gene']
                
                # get actual hit from checked
                try:
                    gene = checked[(checked['GenomeID'] == genome) & (checked['protein'] == hit)].Gene.values[0]
                    if gene == actual_hit:
                        log.write(f'PASS: {hit} is {actual_hit} \n')
                    else:
                        log.write(f'FAIL: expecting {actual_hit} but got {gene} for {hit} \n')
                except:
                    log.write(f'need to update: {hit}, expecting {actual_hit} but got {checked[checked['protein'] == hit].Gene} \n')

            log.write('\n')

def main() -> None:
    test(checked=pd.read_csv('../results/final/nif_final.csv'))

if __name__ == "__main__":
    main()