import pandas as pd
from Bio import SearchIO
import re
import glob
import sys

# get folder
dir = sys.argv[1]

# Get taxonomy from GTDB, https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/
    # NCBI taxonomy from metadata files

GTDB_taxonomy = pd.read_csv('GTDB_taxonomy.tsv', header = None, sep = '\t', names=['GenomeID', 'GTDB'])


'''Parse HMM output and grab hits; hits must meet the following criteria
- positive bit score
- full sequence evalue must be significant (<0.01)
- best domain evalue should be significant (<0.01)
    - otherwise flagged for manual review to check if it is distant homolog or just short repeats

Store hits in feather file'''

# Create empty lists
result_target = []
query_id = []
hit_id = []
evalue = []
best_domain_evalue = []
bitscore = []
bias = []
location = []
orientation = []
alens = []
slength = []
flag1 = []
flag2 = []

def append_hit(genomeID, gene, item):
    result_target.append(genomeID)
    query_id.append(gene)
    hit_id.append(item.id)
    evalue.append(item.evalue)
    best_domain_evalue.append(item.hsps[0].evalue)
    bitscore.append(item.bitscore)
    bias.append(item.bias)
    s = r'# ([0-9]+) # ([0-9]+) # (-?\d)'
    location.append(re.match(s, item.description).group(1) + "-" + re.match(s, item.description).group(2))
    orientation.append(item.description.split('# ')[3])  # default orientation if not specified
    # grab full alignment length (need to sum all domains)
    alen = 0
    for domain in item.hsps:
        alen += domain.aln_span
    alens.append(alen)
    slength.append(int(re.match(s, item.description).group(2))-int(re.match(s, item.description).group(1)))

# Parse through files in output directory

for file in glob.glob(f'../results/{dir}/hmmsearch_results/*.out'):

    # RegEx for the GenomeID (double checking that file is really a genome)
    
    if dir != 'Zostera':
        s = r'([\w]+_[\w]+_[\d]+\.[\d])' # GTDB
    else:
        s = r'([\w]+)_nif' # Zostera
    
    genomeID = re.search(s, file).group()

    # Parse file using SearchIO/HmmerIO
    for result in SearchIO.parse(file, 'hmmer3-text'):
        for item in result.hits:

            # grab gene name
            s = r'([a-zA-Z]+)' # ex. arsA
            gene = re.findall(s, result.id)[0]

            # Check for positive bitscore and append the data to the corresponding lists
            if item.bitscore > 0 and item.evalue < 0.01:
                try:
                    item.hsps[0].evalue
                    append_hit(genomeID, gene, item)
            
                    # check if full seq and best domain e-val are significant
                    if item.hsps[0].evalue < 0.01:
                        flag1.append(0)
                    else:
                        # check if "full sequence Eval is sig but best domain is not, keep only if the target sequence "a multidomain remote homolog; but be wary, and watch out for the case where it’s just a repetitive sequence"
                        flag1.append(1)

                    # check if bitscore >> bias (same order of magnitude) as bitscore
                    if item.bias != 0 and item.bitscore/item.bias > 10:
                        flag2.append(0)
                    else:
                        flag2.append(1)

                except Exception as e:
                    # print(f"Error occurred while processing {genomeID}, {gene}, {item.id}: {e}")
                    pass

# create and store dataframe
hits = pd.DataFrame({'GenomeID': result_target, 'Gene': query_id, 'Hit': hit_id, 
                     'E-value': evalue, 'Best Domain E-value': best_domain_evalue, 'Bit Score': bitscore, 'Bias': bias,
                     'Location': location, 'Orientation': orientation, 'Alignment Length': alens, 'Sequence Length': slength, 
                     'Flag_Eval': flag1, 'Flag_Bias': flag2})

# add taxonomy info
hits = pd.merge(hits, GTDB_taxonomy, on = "GenomeID", how = "left")

# save as feather file
hits.to_feather(f'../results/{dir}/hits.feather')
