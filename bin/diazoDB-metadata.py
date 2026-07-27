# plot operon-org of annotated genes
import requests
import pandas as pd
import glob
import re
import os
import sys
import json
import ast
import argparse
import warnings
from pathlib import Path

from Bio import SeqIO
from pygenomeviz import GenomeViz
from helper import default_proteins_dir

warnings.filterwarnings("ignore", category=FutureWarning)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Organize, export, and plot nif genetic neighborhoods. "
            "Gene annotations are derived from microbeannotator."
        )
    )
    parser.add_argument(
        "--prepare",
        action="store_true",
        help="Create operon FASTA inputs for MicrobeAnnotator.",
    )
    parser.add_argument(
        "--data",
        action='store_true',
        help="Pull operon organization from Microbeannotator output.",
    )
    parser.add_argument(
        "--export",
        action='store_true',
        help="Export operon org data to metadata.json for upload to DiazoDB.",
    )

    parser.add_argument(
        "--plot",
        action='store_true',
        help="Plot operon org data."
    )

    parser.add_argument(
        "--proteins_dir",
        help=(
            "Path to GTDB representative protein directories. "
            "Defaults to DIAZODB_PROTEIN_REPS_DIR, then ../protein_faa_reps_latest, "
            "then ../protein_faa_reps_232."
        ),
        default=default_proteins_dir(),
    )
       
    return parser.parse_args()

def ko2gene(ko):
    gene_abv = {'K00532': 'hydA'} # store ko2gene

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

# run before microbeannotator
def get_operon_fasta(results, proteins_dir):
    # grab fasts file for +/-5 genes around nif operon
    output_dir = Path("../operon-org/input-fastas")
    output_dir.mkdir(parents=True, exist_ok=True)
    proteins_path = Path(proteins_dir)

    for _, cluster in results.iterrows(): # iterate through each genome
        contig = cluster['contig']
        genome = cluster['GenomeID']
        operon = cluster['operon']
        cl = cluster['cluster']
        positions = cluster['pos_num']
        if isinstance(positions, str):
            positions = ast.literal_eval(positions)

        # get positions within +/-5 genes of operon
        start = min(positions) - 5
        end = max(positions) + 5
        acc = [contig + '_' + str(p) for p in range(start, end)] # get acc

        # save subsets as fasta
        file = glob.glob(f"{proteins_path}/*/{genome}_protein.faa")[0]

        records = [record for record in SeqIO.parse(file, "fasta") if record.id in set(acc)]
        output = output_dir / (f"{genome}_{contig}_{operon}_{cl}_operon.fasta")
        SeqIO.write(records, output, "fasta")

# run after mircobeannotator
def get_plot_data():
    # organize microbeannotator results
    nif = pd.read_csv('../results/final/nif_final.csv', index_col=[0,1])
    gene_data = pd.DataFrame(columns = ['genome', 'contig', 'query_id', 'gene', 'ko_number', 'start', 'end', 'orientation'])

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
            annot.loc[row[1].Hit, 'orientation'] = row[1].Orientation
        
        annot.reset_index(inplace=True)
        gene_data = pd.concat([gene_data, annot[['genome', 'contig', 'query_id', 'gene', 'ko_number', 'start', 'end', 'orientation']]])

    gene_data.to_csv('../operon-org/operon-org-plot-data.csv')
    return gene_data

# export metadata.json for displaying hover info on diazoDB phylo tree
def export_metadata(gene_data, operons):
    operons.set_index(['genome', 'contig', 'operon', 'cluster'], inplace=True)

    # known regulon genes (in sort order)
    reg_genes = ['nifA', 'nifL', 'nifR', 'nifI', 'nifI1', 'nifI2', 'glnB', 'glnK', 'draT', 'draG']

    metadata = {}
    for idx, cluster in operons.iterrows():
        genome = idx[0]
        contig = idx[1]
        cl = idx[3]
        taxonomy = cluster['GTDB Taxonomy']
        environments = cluster['Isolation Source']
        regulon = cluster['Regulon']
        organism = cluster['Organism']
        group = cluster['Group No']

        # get operon data for plotting on interactive tree
        genes = []
        regulon = []
        cluster_gene_data = gene_data.loc[gene_data['genome'] == genome]
        for gene in cluster_gene_data.iterrows():
            # for each nif cluster, store surrounding gene info as list
            genes.append({'gene_id': gene[1].query_id,
                        'gene_name': gene[1].gene,
                        'start': gene[1].start,
                        'end': gene[1].end,
                        'direction': gene[1].orientation,
                        'product': gene[1].ko_number})

            # add gene to regulon
            if (gene[1].gene in reg_genes or (re.fullmatch(r'([a-z]{3}[R])', gene[1].gene))):
                    # does the 2nd statement return true if match obj found?
                regulon.append(gene[1].gene)

        # add operon info to metadata for each cluster (operon start, end, and genes in operon)
        operon = {'region_start': cluster_gene_data['start'].min(),
                'region_end': cluster_gene_data['end'].max(),
                'genes': genes}

        metadata[f"{organism} | {genome} | {contig} | {cl}"] = {'oganism': organism, 'genome': genome, 
            'taxonomy': taxonomy, 'group': group, 'environment': environments, 'regulon':regulon, 'operon': operon}
        
    with open('../results/final/metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)


def plot(gene_data): # plot operon organization
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

def main() -> None:
    args = parse_args()

    results = pd.read_csv('../results/final/nif_clusters.csv')
    

    if args.prepare:
        print("Preparing operon FASTA inputs for MicrobeAnnotator", flush=True)
        get_operon_fasta(results, proteins_dir=args.proteins_dir)

    if args.data and not args.prepare:
        print("Pulling operon organization data from MicrobeAnnotator output", flush=True)
        gene_data = get_plot_data()

        print("Exporting operon organization to metadata.json", flush=True)
        export_metadata(gene_data, results)

        if args.plot:
            plot(gene_data)

    else:
        gene_data = pd.read_csv('../operon-org/operon-org-plot-data.csv', index_col=[0,1])

        if args.export:
            print("Exporting operon organization to metadata.json", flush=True)
            export_metadata(gene_data, results)

        if args.plot:
            print("Plotting operon organization data from existing CSV", flush=True)
            plot(gene_data)


if __name__ == "__main__":
    main()