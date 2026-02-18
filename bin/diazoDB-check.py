import pandas as pd


# grab file nif hits
nif_hits = pd.read_csv(f'../results/nif_final.csv')

# read test file
test_file = pd.read_csv('diazoDB-checks.csv')
test_file.set_index(['genomeID', 'taxonomy'], inplace=True)

# create log file
with open('diazoDB-check-log.txt', 'w') as log:
    # iterate through genomes in test file
    for genome in test_file.index.get_level_values(0).unique():
        taxonomy = test_file.loc[genome].index.get_level_values(0)[0]
        # write genome, taxonomy to log
        log.write(f'Checking {genome} {taxonomy}\n')

        # compare if correct # genes found
        nif_hits_count = nif_hits[nif_hits['GenomeID'] == genome].shape[0]
        test_file_count = test_file.loc[genome].shape[0]

        if nif_hits_count == test_file_count:
            log.write('PASS: correct number of genes found\n')
        elif nif_hits_count < test_file_count:
            log.write(f'FAIL: missing {test_file_count - nif_hits_count} hits \n')
        elif nif_hits_count > test_file_count:
                extra = [hit for hit in nif_hits[nif_hits['GenomeID'] == genome].Hit.to_list() if hit not in test_file.loc[genome].hit.to_list()]
                log.write(f'FAIL: {nif_hits_count - test_file_count} extra hits: {extra} \n')

        # compare gene annotations
        for index,row in test_file.loc[genome].iterrows():

            # get test values
            hit = row['hit']
            actual_hit = row['gene']
            
            # get actual hit from nif_hits
            try:
                gene = nif_hits[(nif_hits['GenomeID'] == genome) & (nif_hits['Hit'] == hit)].Gene.values[0]
                if gene == actual_hit:
                    log.write(f'PASS: {hit} is {actual_hit} \n')
                else:
                    log.write(f'FAIL: expecting {actual_hit} but got {gene} for {hit} \n')
            except:
                log.write(f'need to update: {hit}, expecting {actual_hit} but got {gene} \n')

        log.write('\n')