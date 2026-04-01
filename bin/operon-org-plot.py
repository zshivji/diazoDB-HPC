# plot operon-org of annotated genes
import requests
import pandas as pd
import glob
import re
import sys
import json

from Bio import SeqIO
from pygenomeviz import GenomeViz

# organize microbeannotator results
nif = pd.read_csv('../results/final/nif_final.csv', index_col=[0,1])

gene_abv = {'K00532': 'hydA'} # store ko2gene
gene_data = pd.DataFrame(columns = ['genome', 'contig', 'query_id', 'gene', 'ko_number', 'start', 'end', 'orientation'])

def ko2gene(ko):
    try:
        return gene_abv[ko]
    except:
        url = f"https://rest.kegg.jp/get/ko:{ko}"
        r = requests.get(url)
        for line in r.text.split("\n"):
            if line.startswith("SYMBOL"):
                gene = line.split()[-1]
                gene_abv[ko] = gene
                return gene

# grab annotation files
for file in glob.glob(f"../operon-org/microbeannotator/annotation_results/*.annot"):
    # convet ko_number to gene abv
    annot = pd.read_csv(file, sep = '\t')
    annot['gene'] = annot['ko_number'].apply(ko2gene)
    
    # grab genome and contig from file name
    str = r'([A-Z_]+[A-Z0-9.]+)_([A-Z_]*[A-Z0-9.]+)'
    genome = re.search(str, file).group(1)
    annot['genome'] = genome
    contig = re.search(str, file).group(2)
    annot['contig'] = contig

    # grab start, end, orientation of each gene
    input = glob.glob(f"../operon-org/input-fastas/{genome}_{contig}_operon.fasta")[0]
    for result in SeqIO.parse(input, "fasta"):
        start = int(result.description.split('# ')[1])
        end = int(result.description.split('# ')[2])
        orientation = int(result.description.split('# ')[3])
        if result.id in annot.index:
            annot.loc[result.id, 'start'] = start
            annot.loc[result.id, 'end'] = end
            annot.loc[result.id, 'orientation'] = orientation
    
    for row in nif.loc[(genome, contig)].iterrows():
        annot.loc[row[1].Hit, 'gene'] = row[1].Gene
        annot.loc[row[1].Hit, 'start'] = int(row[1].Location.split('-')[0])
        annot.loc[row[1].Hit, 'end'] = int(row[1].Location.split('-')[1])
        annot.loc[row[1].Hit, 'orientation'] = 1 # BUG
    
    annot.reset_index(inplace=True)
    gene_data = pd.concat([gene_data, annot[['genome', 'contig', 'query_id', 'gene', 'ko_number', 'start', 'end', 'orientation']]])

gene_data.to_csv('../operon-org/operon-org-plot-data.csv')

# update metadata.json to store gene data
with open('../results/final/metadata.json', 'a') as file:
    metadata = json.load(file)

    for cluster in metadata.items():
        print(cluster)
        cluster_gene_data = gene_data.loc[gene_data['genome'] == cluster[1]['genome']]
        print(cluster_gene_data)

        genes = []
        for gene in cluster_gene_data.iterrows():
            # for each nif cluster, store surrounding gene info as list
            genes.append({'gene_id': gene[1].query_id,
                          'gene_name': gene[1].gene,
                          'start': gene[1].start,
                          'end': gene[1].end,
                          'direction': gene[1].orientation,
                          'product': gene[1].ko_number})
        # add operon info to metadata for each cluster (operon start, end, and genes in operon)
        operon = {'region_start': cluster_gene_data['start'].min(),
                  'region_end': cluster_gene_data['end'].max(),
                  'genes': genes}
        
        metadata[cluster[1]['operon']] = operon

    # save updated metadata
    json.dump(metadata, file)

# plot operon organization
to_plot = sys.argv[0]

if to_plot:
    gene_data = pd.read_csv('../operon-org/operon-org-plot-data.csv', index_col=[0,1])
    genome_list = gene_data.index.get_level_values(0).unique().tolist()

    gv = GenomeViz()

    for genome in genome_list:
        genes = gene_data.loc[(genome)].gene.to_list()
        starts = gene_data.loc[(genome)].start.to_list()
        ends = gene_data.loc[(genome)].end.to_list()
        orientations = gene_data.loc[(genome)].orientation.to_list()

        track = gv.add_feature_track(genome, (int(min(starts)), int(max(ends))))
        for idx, gene in enumerate(genes):
            if gene == 'nifH':
                color = 'blue'
            elif gene == 'nifD':
                color = 'red'
            elif gene == 'nifK':
                color = 'green'
            elif gene == 'nifB':
                color = 'purple'
            elif gene == 'nifE':
                color = 'orange'
            elif gene == 'nifN':
                color = 'pink'
            else:
                color = 'grey'
            
            track.add_feature(
                int(starts[idx]),
                int(ends[idx]),
                int(orientations[idx]),
                plotstyle='bigarrow',
                fc = color,
                lw = 1,
                label = gene,
                text_kws=dict(rotation=0, vpos="center", hpos="center"))

    gv.savefig("../operon-org/ALL.png")
