import pandas as pd
from Bio import SeqIO
import glob

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from cluster_pos import cluster_pos

# make sure gene clusters have at least nifHDK
def gene_check(genes):
    if genes.__contains__('nifH'):
        if genes.__contains__('nifD'):
            if genes.__contains__('nifK'):
                return True

# multi-index to cluster by genome, contig
nif = pd.read_csv(f'../results/tmp/nif_rescheck_nofilt.csv')
nif.reset_index(inplace = True)
nif.set_index(['GenomeID', 'contig'], inplace = True)
nif.sort_index(inplace = True)
nif.drop_duplicates(inplace = True)

# filter for genome, contig with at least 3 unique genes (nifHDKENB)
filtered_nif = nif.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Gene'].nunique() >= 3)

# make sure these 3 unique genes are not the same hit (i.e. not the same gene in reference genome)
filtered_nif2 = filtered_nif.groupby(level=['GenomeID', 'contig']).filter(lambda x: x['Hit'].nunique() >= 3)

genomes_to_keep = []
# iterate through each genome and contig
for genome in filtered_nif2.index.get_level_values(0).unique(): # iterate through each genome
    for contig in filtered_nif2.loc[genome].index.get_level_values(0).unique(): # iterate through each contig

        tmp = filtered_nif2.loc[(genome, contig)]

        # only keep numbers that have clusters >= 3
        pos_clusters = cluster_pos(tmp.Hit.unique(), 20)

        # for each cluster, find the best combination of genes (min e-value)
        for ind, cl in enumerate(pos_clusters):
            pos = [contig + '_' + str(p) for p in cl]
            no_pos = len(pos)
            
            # need at least 3 genes to continue
            if no_pos < 3:
                continue

            # # only keep hits that are in the cluster --> doesn't work for nifB
            # tmp2 = tmp[tmp.Hit.isin(pos)]

            # # only keep hits that are in the cluster
            # tmp2 = tmp[tmp.Hit.isin(pos)].reset_index()

            # check if all genes are present
            if gene_check(tmp.residue_match.to_list()):
                # get index
                items = [(genome, contig, hit) for hit in tmp.Hit]
                genomes_to_keep.extend(items)

# filter for genomes to keep
filtered_nif2.set_index(['Hit'], append = True, inplace = True)
filtered_nif2 = filtered_nif2.loc[genomes_to_keep]
filtered_nif2.sort_index(inplace = True)

#clean up cols
filtered_nif2['Gene'] = filtered_nif2['residue_match']
filtered_nif2 = filtered_nif2[['Gene', 'E-value', 'Bit Score', 'Location', 'Orientation', 'Alignment Length', 'Sequence Length', 'GTDB']]
filtered_nif2.drop_duplicates(inplace = True)

# export csv with each gene as individual rowfiltered_nif2.to_feather(f'../results/final/nif_final.feather')
filtered_nif2.to_csv(f'../results/final/nif_final.csv')

# export csv grouped by genome
    # works but can simplify code once trees are built w hit not GenomeID

def get_hit(rec):
    return rec.description.split(' ')[-1]

def get_cluster(file):
    clusters = {}
    clusters_fasta = list(SeqIO.parse(file, 'fasta'))

    index = 0
    while index < len(clusters_fasta)-1:
        rec = clusters_fasta[index]
        if rec.seq == '': # cluster header
            cluster = get_hit(clusters_fasta[index+1]) # hit
            clause = True
            members = []
            i = 1
            while clause & (index+i < len(clusters_fasta)): # add members (until next cluster header)
                if clusters_fasta[index+i].seq != '':
                    members.append(get_hit(clusters_fasta[index+i]))
                    i+=1
                elif clusters_fasta[index+i].seq == '':
                    clause = False
            index += i

            clusters[cluster] = members

    return clusters

nif = pd.read_csv('../results/final/nif_final.csv')
nif.reset_index(inplace = True)
nif.set_index(['Hit'], inplace=True)

# add nif groups (based on nifH--for now)

nif['Group'] = ''

#genes = {'H': 'H', 'D': 'D_noOut'}
genes = {'H': 'H'}

for gene, file in genes.items():
    # get clustered datapoints
    clusters = get_cluster(f'trees/nif{file}/clustered_nif{file}_all_seqs.fasta') 

    # assign group 
    for group in ['1', '2', '3', '4a', '4c', '3anfvnf']:
        lines = []
        hits = []
        with open(f'bin/nif_groups/nif{gene}_group{group}.txt','r') as f:
            lines = f.read().splitlines()
            for line in lines:
                hit = '_'.join(line.split('|')[-1].split(' '))[:-1]
                hits.append(hit) # reformat "hits" to match nif index
                hits.extend(clusters[hit]) # add clustered hits to list of hits to update
        for hit in hits:
            nif.loc[nif.index == hit, 'Group'] = f'Group {group}'

# load GTDB metadata
GTDB_metadata = pd.read_csv('GTDB_metadata.gz', sep = '\t', usecols=['accession', 'gtdb_taxonomy', 'ncbi_taxonomy', 
                                                'ncbi_taxonomy_unfiltered', 'ncbi_country', 'ncbi_isolation_source'])

# def for sorting gene order
def sort_order(g):
    order = ['H', 'D', 'K', 'B', 'E', 'N']
    return order.index(g)

# get gene set (i.e. HDKEN)
nif['Gene'] = nif['Gene'].str.replace('nif', '')
nif['Gene set'] = nif[['GenomeID','contig','Gene']].groupby(['GenomeID','contig'])['Gene'].transform(lambda x: ''.join(sorted(x, key=sort_order)))
# get contig positions (i.e. 1,2,3)
nif['index'] = nif['Hit'].str.split('_').str[-1]
nif['Position'] = nif[['GenomeID','contig','index']].groupby(['GenomeID','contig'])['index'].transform(lambda x: ','.join(sorted(x)))

# group by genome, contig and save
nif = nif[['GenomeID', 'contig', 'Gene set', 'Position', 'GTDB', 'Location' ,'Orientation', 'Group']]
nif = nif.groupby(['GenomeID','contig']).first()
nif.reset_index(level=['GenomeID','contig'], inplace = True)

#merge dataframes to include accession, metadata, and sequences (for filtering)
nif = pd.merge(nif, GTDB_metadata, left_on = 'GenomeID', right_on = 'accession', how = 'left').drop(columns = ['accession', 'GTDB'])

# add summary genomic location (contig:min-max of cluster, like NCBI format)
nif['Organism'] = nif['gtdb_taxonomy'].str.split(';').str[-1].str.split('__').str[-1]
nif['Regulon'] = ''
nif['PredGrowthTemp'] = ''
nif.rename(columns = {'GenomeID': 'Genome', 'contig': 'Contig', 'Gene set': 'Nitrogenase Set', 
                      'Group':'Group No', 'gtdb_taxonomy': 'GTDB Taxonomy', 'ncbi_taxonomy': 'NCBI Taxonomy',
                      'ncbi_isolation_source':'Isolation Source'}, inplace = True)
nif.sort_index(inplace = True)

nif.to_csv(f'../results/final/nif_genomes.csv')

# export metadata.json for displaying hover info on diazoDB phylo tree
metadata = {}

for cluster in nif.iterrows():
    genome = cluster[1]['Genome']
    taxonomy = cluster[1]['GTDB Taxonomy']
    environments = cluster[1]['Isolation Source']
    regulon = cluster[1]['Regulon']
    operon = ''

    metadata[cluster[1]['Organism']] = {'genome': genome, 'taxonomy': taxonomy,
                                        'environments': environments, 'regulon':regulon, 'operon': operon}
    
import json
with open('../results/final/metadata.json', 'w') as f:
    json.dump(metadata, f)

# export fasta files (use .copy() to avoid SettingWithCopyWarning)
nif = pd.read_feather(f'../results/final/nif_final.feather')

nifH = nif[(nif.Gene == 'nifH')].copy()
nifD = nif[(nif.Gene == 'nifD')].copy()
nifK = nif[(nif.Gene == 'nifK')].copy()
nifB = nif[(nif.Gene == 'nifB')].copy()
nifE = nif[(nif.Gene == 'nifE')].copy()
nifN = nif[(nif.Gene == 'nifN')].copy()
nifNB = nif[(nif.Gene == 'nifNB')].copy()

# anfH = nif[(nif.Gene == 'anfH')].copy()
# anfD = nif[(nif.Gene == 'anfD')].copy()
# anfK = nif[(nif.Gene == 'anfK')].copy()
# anfO = nif[(nif.Gene == 'anfO')].copy()
# anfG = nif[(nif.Gene == 'anfG')].copy()

# vnfH = nif[(nif.Gene == 'vnfH')].copy()
# vnfD = nif[(nif.Gene == 'vnfD')].copy()
# vnfK = nif[(nif.Gene == 'vnfK')].copy()
# vnfO = nif[(nif.Gene == 'vnfO')].copy()
# vnfG = nif[(nif.Gene == 'vnfG')].copy()


gene_list = [nifH, nifD, nifK, nifB, nifE, nifN, nifNB]
gene_names = ['nifH', 'nifD', 'nifK', 'nifB', 'nifE', 'nifN', 'nifNB']

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
    with open("../results/final/fastas/final_" + name + ".fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

print("done", flush=True)