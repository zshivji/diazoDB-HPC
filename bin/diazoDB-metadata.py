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
    parser.add_argument(
        "--clusters_file",
        type=Path,
        default=Path("../results/final/nif_clusters.csv"),
        help="Nif cluster table used to select operon neighborhoods.",
    )
    parser.add_argument(
        "--nif_final_file",
        type=Path,
        default=Path("../results/final/nif_final.csv"),
        help="Detailed nif calls used to restore DiazoDB gene assignments.",
    )
    parser.add_argument(
        "--operon_dir",
        type=Path,
        default=Path("../operon-org"),
        help="Job-local operon working directory.",
    )
    parser.add_argument(
        "--metadata_file",
        type=Path,
        default=Path("../results/final/metadata.json"),
        help="Destination for exported operon metadata JSON.",
    )
    parser.add_argument(
        "--plot_file",
        type=Path,
        default=Path("../operon-org/ALL.png"),
        help="Destination for the operon organization plot.",
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
def get_operon_fasta(results, proteins_dir, operon_dir):
    # grab fasts file for +/-5 genes around nif operon
    output_dir = Path(operon_dir) / "input-fastas"
    output_dir.mkdir(parents=True, exist_ok=True)
    proteins_path = Path(proteins_dir)

    for _, cluster in results.iterrows(): # iterate through each genome
        contig = cluster['contig']
        genome = cluster['GenomeID']
        operon = cluster['operon']
        cl = cluster['cluster']
        # Database-build tables contain ``pos_num``; external DiazoDB tables
        # contain the equivalent ``proteins`` list (contig_gene_number).
        positions = cluster.get('pos_num')
        if positions is None or (isinstance(positions, float) and math.isnan(positions)):
            proteins = cluster.get('proteins', [])
            if isinstance(proteins, str):
                proteins = ast.literal_eval(proteins)
            positions = [int(str(protein).rsplit('_', 1)[1]) for protein in proteins]
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
def get_plot_data(nif_final_file, clusters_file, operon_dir):
    # organize microbeannotator results
    nif = pd.read_csv(nif_final_file)
    if 'operon' not in nif.columns:
        clusters = pd.read_csv(clusters_file)
        nif = nif.merge(
            clusters[['GenomeID', 'contig', 'cluster', 'operon']],
            on=['GenomeID', 'contig', 'cluster'],
            how='left',
            validate='many_to_one',
        )
    annots =[]

    # for each nif cluster, store surrounding operon data
    for (genome, contig, operon, cluster), subset in nif.groupby(['GenomeID', 'contig', 'operon', 'cluster'], sort=False):

        file = Path(operon_dir) / "microbeannotator" / "annotation_results" / f"{genome}_{contig}_{operon}_{cluster}_operon.fasta.annot"
        annot = pd.read_csv(file, sep = '\t', index_col = 'query_id')

        # convert ko_number to gene abv
        annot['ko_number'] = annot['ko_number'].apply(ko2gene)
        # update gene column to use ko_number if available, otherwise use protein_id
        annot['gene'] = annot['ko_number'].combine_first(annot['protein_id'])

        # for each annotated gene, grab start, end, and orientation from fasta header
        fasta_file = Path(operon_dir) / "input-fastas" / f"{genome}_{contig}_{operon}_{cluster}_operon.fasta"
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
    gene_data = pd.concat(annots, ignore_index=True) if annots else pd.DataFrame(
        columns=['query_id', 'genome', 'contig', 'operon', 'cluster',
                 'gene', 'start', 'end', 'orientation']
    )
    gene_data.to_csv(Path(operon_dir) / 'operon-org-plot-data.csv', index=False)
    return gene_data

# need tree built with final data to get group info
    # appends Grou info onto exisiting nif_final.csv and nif_clusters.csv
def get_group():
    nif = pd.read_csv('../results/final/nif_final.csv')

    # start with assigning nifH clusters
    gene = 'H'

    # get clustered datapoints
    tree_clusters = pd.read_csv(f'../trees/nif{gene}/nif{gene}_anf{gene}_vnf{gene}_clustered.fasta.tsv', 
                                sep = '\t', 
                                header = None,
                                names = ['rep', 'acc']) 

    # assign group 
    for group in ['1', '2', '3', '4a', '4c', '3anfvnf']:
        lines = []
        hits = []
        # for each group, find all matching hits in nif_final.csv
        with open(f'nif_groups/nif{gene}_group{group}.txt','r') as f:
            lines = f.read().splitlines()
            for line in lines:
                hit = '_'.join(line.split('|')[-1].strip().replace("'", "").split(' '))
                hits.append(hit) # reformat "hits" to match nif index
                hits.extend(tree_clusters.loc[tree_clusters['rep'] == hit, 'acc'].to_list()) # add clustered hits to list of hits to update
                if group == '3anfvnf':
                    hits.extend(nif.loc[nif['Gene'].isin(['anfH', 'vnfH']), 'protein'].to_list()) # add anfH/vnfH to group 3anfvnf
        
        # Apply the nifH/anfH/vnfH group to every gene in each matched cluster
        cluster_cols = ['GenomeID', 'contig', 'cluster', 'operon']

        grouped_clusters = nif.loc[nif['protein'].isin(hits), cluster_cols].drop_duplicates()
        cluster_index = pd.MultiIndex.from_frame(grouped_clusters)
        nif_index = pd.MultiIndex.from_frame(nif[cluster_cols])
        nif.loc[nif_index.isin(cluster_index), 'Group'] = f'Group {group}'
    
    # export updated nif_final.csv with group info
    nif.to_csv('../results/final/nif_final.csv', index=False)

    # export updated nif_clusters.csv with group info
    clusters = pd.read_csv('../results/final/nif_clusters.csv')
    clusters = clusters.drop(columns=['Group', 'Group No', 'Group_x', 'Group_y'], errors='ignore') # make sure no duplicate Group columns exist before merging
    # add Group col to nif_clusters.csv by matching rows GenomID, contig, cluster, and operon to nif_final.csv
        # how='left' --> keep all rows in nif_clusters.csv, even if no match in nif_final.csv
        # validate='many_to_one' --> each row in nif_clusters.csv should match at most one row in nif_final.csv
    clusters = clusters.merge(nif[['GenomeID', 'contig', 'cluster', 'operon', 'Group']].drop_duplicates(), 
                              on=['GenomeID', 'contig', 'cluster', 'operon'], how='left', validate='many_to_one')
    clusters.to_csv('../results/final/nif_clusters.csv', index=False)

# export metadata.json for displaying hover info on diazoDB phylo tree
def export_metadata(gene_data, operons, metadata_file):
    indexed_gene_data = gene_data.set_index(
        ['genome', 'contig', 'operon', 'cluster']
    )

    # known regulon genes (in sort order)
    reg_genes = ['nifA', 'nifL', 'nifR', 'nifI', 'nifI1', 'nifI2', 'glnB', 'glnK', 'draT', 'draG', 'cnfR']

    metadata = {}
    for _, cluster in operons.iterrows():
        genome = cluster['GenomeID']
        contig = cluster['contig']
        cl = cluster['cluster']
        operonID = cluster['operon']
        taxonomy = cluster.get('GTDB Taxonomy', cluster.get('GTDB', ''))
        environments = cluster.get('Isolation Source', '')
        organism = cluster.get('Organism', genome)
        group = cluster.get('Group No', cluster.get('Group', ''))

        # get operon data for plotting on interactive tree
        genes = []
        regulon = []
        cluster_gene_data = indexed_gene_data.loc[
            [(genome, contig, operonID, cl)]
        ]
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
        
    Path(metadata_file).parent.mkdir(parents=True, exist_ok=True)
    with open(metadata_file, 'w') as f: # overwrites existing metadata.json
        json.dump(json_safe(metadata), f, indent=2, allow_nan=False)


def plot(gene_data, plot_file): # plot operon organization
    gv = GenomeViz()

    for genome, genome_data in gene_data.groupby('genome', sort=False):
        genes = genome_data.gene.to_list()
        starts = genome_data.start.to_list()
        ends = genome_data.end.to_list()
        orientations = genome_data.orientation.to_list()

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

    Path(plot_file).parent.mkdir(parents=True, exist_ok=True)
    gv.savefig(plot_file)

def main() -> None:
    args = parse_args()

    if args.prepare:
        results = pd.read_csv(args.clusters_file)
        print("Preparing operon FASTA inputs for MicrobeAnnotator", flush=True)
        get_operon_fasta(
            results,
            proteins_dir=args.proteins_dir,
            operon_dir=args.operon_dir,
        )

    if args.data:
        print("Pulling operon organization data from MicrobeAnnotator output", flush=True)
        gene_data = get_plot_data(
            nif_final_file=args.nif_final_file,
            clusters_file=args.clusters_file,
            operon_dir=args.operon_dir,
        )

        get_group() # make sure nif group info is appended to nif_final.csv and nif_clusters.csv
        results = pd.read_csv(args.clusters_file)

        print("Exporting operon organization to metadata.json", flush=True)
        export_metadata(gene_data, results, metadata_file=args.metadata_file)

        if args.plot:
            plot(gene_data, plot_file=args.plot_file)

    elif args.export or args.plot:
        get_group() # make sure nif group info is appended to nif_final.csv and nif_clusters.csv
        results = pd.read_csv(args.clusters_file)
        gene_data = pd.read_csv(args.operon_dir / 'operon-org-plot-data.csv')

        if args.export:
            print("Exporting operon organization to metadata.json", flush=True)
            export_metadata(gene_data, results, metadata_file=args.metadata_file)

        if args.plot:
            print("Plotting operon organization data from existing CSV", flush=True)
            plot(gene_data, plot_file=args.plot_file)


if __name__ == "__main__":
    main()
