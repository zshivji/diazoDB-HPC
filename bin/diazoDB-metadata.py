# plot operon-org of annotated genes
import requests
import pandas as pd
import glob
from functools import lru_cache
import re
import os
import sys
import json
import ast
import argparse
import warnings
import math
import numbers
from pathlib import Path

from Bio import SeqIO
from pygenomeviz import GenomeViz
from helper import default_proteins_dir

warnings.filterwarnings("ignore", category=FutureWarning)


def json_safe(value):
    """Convert values unsupported by strict JSON to JSON-compatible values.
    This prevents errors when js renders metadata.json on diazoDB home page."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, numbers.Real) and math.isnan(value):
        return None
    return value

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

KO_OVERRIDES = {'K00532': 'hydA'} # store ko2gene

@lru_cache(maxsize=None)
def ko2gene(ko):
    # if known override exists, return that instead of querying KEGG
    if ko in KO_OVERRIDES:
        return KO_OVERRIDES[ko]

    url = f"https://rest.kegg.jp/get/ko:{ko}"
    r = requests.get(url)
    for line in r.text.splitlines():
        if line.startswith("SYMBOL"):
            return line.split()[-1]

    return ko # fallback

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
    nif = pd.read_csv('../results/final/nif_final.csv')
    annots =[]

    # for each nif cluster, store surrounding operon data
    for (genome, contig, operon, cluster), subset in nif.groupby(['GenomeID', 'contig', 'operon', 'cluster'], sort=False):

        file = f"../operon-org/microbeannotator/annotation_results/{genome}_{contig}_{operon}_{cluster}_operon.fasta.annot"
        annot = pd.read_csv(file, sep = '\t', index_col = 'query_id')

        # convert ko_number to gene abv
        annot['ko_number'] = annot['ko_number'].apply(ko2gene)
        # update gene column to use ko_number if available, otherwise use protein_id
        annot['gene'] = annot['ko_number'].combine_first(annot['protein_id'])

        # for each annotated gene, grab start, end, and orientation from fasta header
        fasta_file = f"../operon-org/input-fastas/{genome}_{contig}_{operon}_{cluster}_operon.fasta"
        fasta = SeqIO.parse(fasta_file, "fasta")

        for seq in fasta:
            if seq.id not in annot.index:
                continue
            else:
                annot.loc[seq.id, 'start'] = int(seq.description.split('# ')[1])
                annot.loc[seq.id, 'end'] = int(seq.description.split('# ')[2])
                annot.loc[seq.id, 'orientation'] = int(seq.description.split('# ')[3])

        # for diazoDB annotated nif genes, replace with diazoDB results
            # accounts for incorrect annotations from microbeannotator (e.g. anfH annotated as nifH)
        for _, row in subset.iterrows():
            annot.loc[row.protein, 'gene'] = row.Gene
            annot.loc[row.protein, 'start'] = int(row.Location.split('-')[0])
            annot.loc[row.protein, 'end'] = int(row.Location.split('-')[1])
            annot.loc[row.protein, 'orientation'] = row.Orientation
            if row.Gene == 'anfO':
                print(row.protein)

        # add info at end, after all rows added (microbeannotator + diazoDB -- microbeannotator may have missed some diazoDB hits (anfO, nifB, etc.))
        annot['genome'] = genome
        annot['contig'] = contig
        annot['operon'] = operon
        annot['cluster'] = cluster

        annot = annot[['genome', 'contig', 'operon', 'cluster', 'gene', 'start', 'end', 'orientation']]
        annot.reset_index(inplace=True)
        annots.append(annot)

    #gene_data = pd.DataFrame(columns = ['genome', 'contig', 'query_id', 'gene', 'ko_number', 'start', 'end', 'orientation'])
    gene_data = pd.concat(annots, ignore_index=True)
    gene_data.to_csv('../operon-org/operon-org-plot-data.csv', index=False)
    return gene_data

# export metadata.json for displaying hover info on diazoDB phylo tree
def export_metadata(gene_data, operons):

    # known regulon genes (in sort order)
    reg_genes = ['nifA', 'nifL', 'nifR', 'nifI', 'nifI1', 'nifI2', 'glnB', 'glnK', 'draT', 'draG']

    metadata = {}
    for _, cluster in operons.iterrows():
        genome = cluster['GenomeID']
        contig = cluster['contig']
        cl = cluster['cluster']
        operonID = cluster['operon']
        taxonomy = cluster['GTDB Taxonomy']
        environments = cluster['Isolation Source']
        regulon = cluster['Regulon']
        organism = cluster['Organism']
        group = cluster['Group No']

        # get operon data for plotting on interactive tree
        genes = []
        regulon = []
        cluster_gene_data = gene_data.loc[gene_data['genome'] == genome]
        for _, gene in cluster_gene_data.iterrows():
            # for each nif cluster, store surrounding gene info as list
            genes.append({'gene_id': gene.query_id,
                        'gene_name': gene.gene,
                        'start': gene.start,
                        'end': gene.end,
                        'direction': gene.orientation})

            # add gene to regulon
            try:
                if re.fullmatch(r'([a-z]{3}R)', gene.gene):
                    regulon.append(gene.gene)
            except:
                pass
            if (gene.gene in reg_genes):
                regulon.append(gene.gene)

        # add operon info to metadata for each cluster (operon start, end, and genes in operon)
        operon = {'region_start': cluster_gene_data['start'].min(),
                'region_end': cluster_gene_data['end'].max(),
                'genes': genes}

        metadata[f"{organism} | {cl} | {genome} | {contig} | {operonID}"] = {'organism': organism, 'genome': genome, 
            'taxonomy': taxonomy, 'group': group, 'environment': environments, 'regulon':regulon, 'operon': operon}
        
    with open('../results/final/metadata.json', 'w') as f:
        json.dump(json_safe(metadata), f, indent=2, allow_nan=False)


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
        gene_data = pd.read_csv('../operon-org/operon-org-plot-data.csv')

        if args.export:
            print("Exporting operon organization to metadata.json", flush=True)
            export_metadata(gene_data, results)

        if args.plot:
            print("Plotting operon organization data from existing CSV", flush=True)
            plot(gene_data)


if __name__ == "__main__":
    main()
