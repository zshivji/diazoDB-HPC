import pandas as pd
import numpy as np
import glob
import json
import os
import subprocess
from tqdm import tqdm
import argparse

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from helper import load_nif, filter_groups_by_unique_counts, get_seq, default_proteins_dir
from test import test

# get args (reload_fasta)
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make and align temporary fasta files based on filtered hits."
            "Check for conserved residues in alignments."
            "Update nif gene assignments based on conserved residue matching."
        )
    )
    parser.add_argument(
        "--reload_fasta",
        help="Whether or not to reload the fasta files (slow but necessary if hmmsearch/parse_hmmp.py changed).",
        action='store_true'
    )

    parser.add_argument(
        "--align_fasta",
        help="Whether or not to align the fasta files (slow but necessary if hmmsearch/parse_hmmp.py changed). Automatically defaults to True if --reload_fasta is specified.",
        action='store_true'
    )

    parser.add_argument(
        "--config_file",
        help="Path to the JSON configuration file containing reference sequences and important residues. Default bin/nif-config.json",
        default='nif-config.json'
    )

    parser.add_argument(
        "--proteins_dir",
        help=(
            "Root containing <group>/<genome>_protein.faa files. "
            "The API runner supplies its uploaded proteins and packaged references."
        ),
        default=default_proteins_dir(),
    )

    parser.add_argument(
        "--min_genes",
        type=int,
        default=3,
        help="Minimum number of nif genes required to be considered a tophit cluster (default: 3)."
    )

    parser.add_argument(
        "--test",
        action='store_true',
        help="Whether or not to run the final results against the test file."
    )
    parser.add_argument(
        "--hits_file",
        help=(
            "Optional parsed HMM hits CSV. Defaults to the repository-wide "
            "archaea_hits.csv and bacteria_hits.csv files."
        ),
    )
    parser.add_argument(
        "--results_dir",
        default="../results/tmp",
        help="Directory for temporary FASTA, alignment, and residue-check files.",
    )
    parser.add_argument(
        "--final_dir",
        default="../results/final",
        help="Directory for nif_final.csv, nif_clusters.csv, and final FASTAs.",
    )
    parser.add_argument(
        "--ref_seq_ids",
        default="../results/ref_seq.ids",
        help="seqtk ID file identifying the conserved-residue reference genome.",
    )
    parser.add_argument(
        "--skip_metadata",
        action="store_true",
        help="Do not read or join repository taxonomy metadata (used by API runner jobs).",
    )

    parser.add_argument(
        "--external",
        action="store_true",
        help="Export a cleaned-up nif_final.csv and nif_clusters.csv for external DiazoDB users.",
    )

    return parser.parse_args()

# Note: runtime data (nif table, config, per-gene frames) are loaded in `main`.
# Helper functions below accept explicit inputs and return results.

def write_fasta(
    gene_list: dict,
    config: dict,
    proteins_dir: str = "../protein_faa_reps_232",
    results_dir: str = "../results/tmp",
):
    """Get fasta sequences for each gene & export to fasta of domains.
    """
    os.makedirs(results_dir, exist_ok=True)
    for name, gene in gene_list.items():
        print('getting ' + name + ' fasta', flush=True)
        # filter out hits that don't meet length criteria (domain is too short, pre-filtering decreases runtime)
        minl = config[name]['min-length']
        gene = gene[abs(gene['Domain Start'] - gene['Domain End']) >= minl]
        gene_list[name] = gene

        records = []
        for genome, hit in gene.iterrows():
            seq = get_seq(
                genome,
                hit.Hit[:-2],
                description=hit.Hit,
                start=hit['Domain Start'],
                end=hit['Domain End'],
                dir=proteins_dir,
            )
            records.append(seq)

        print(len(records), "records found for " + name, flush=True)
        with open(os.path.join(results_dir, name + ".fasta"), "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

    gene_list = pd.concat(gene_list.values())
    gene_list.to_csv(os.path.join(results_dir, 'hits_seq.csv'))

    return gene_list

def aln_fasta(
    genes: dict,
    config: dict,
    results_dir: str = "../results/tmp",
):
    # align fasta files
    splits_dir = os.path.join(results_dir, 'fasta_splits')
    os.makedirs(splits_dir, exist_ok=True)
    for name in genes['Gene'].unique():
        print("aligning "+ name, flush=True)
        tmp = genes[genes['Gene']==name].copy()
        num = int(tmp.shape[0]/200)+1 # how many splits
        fasta_path = os.path.join(results_dir, f"{name}.fasta")
        split_prefix = os.path.join(splits_dir, f"{name}_split")
        ref_fasta_path = f"/resnick/groups/enviromics/zahra/diazoDB-HPC/results/reference/{name}.fasta" # placeholder for ref fasta path (for diazoDB and API runner, will always point to AV)
        # split fasta file into smaller files for alignment
        subprocess.run(
            ["seqtk", "split", "-n", str(num), split_prefix, fasta_path],
            check=True,
        )
        # for each split...
        for i in range(num):
            split_path = f"{split_prefix}.{str(i+1).zfill(5)}.fa"
            alignment_path = f"{split_prefix}.{str(i+1).zfill(5)}.aln"
            
            # add reference sequences to split file
            ref = SeqRecord(Seq(config[name]['ref_seq']), id='reference', description=name)
            with open(split_path, "a") as f:
                SeqIO.write(ref, f, "fasta")

            # align split file
            with open(alignment_path, "wb") as alignment_file:
                subprocess.run(
                    ["mafft", "--auto", "--quiet", "--thread", "4", split_path],
                    stdout=alignment_file,
                    check=True,
                )


def check_gene(file, ref_seq, important_residues, residue_scores, passing_score, p=False):
    if len(important_residues) != len(residue_scores):
        raise ValueError(
            f"{file[-20:-16]}: residues and residue_scores must have the same length "
            f"({len(important_residues)} != {len(residue_scores)})"
        )

    residue_to_score = dict(zip(important_residues, residue_scores))
    alignment = AlignIO.read(f"{file}", "fasta")

    for record in alignment:
        if record.id == 'reference':
            aln = record.seq
            break

    # map important cols in aln to ref_seq
    residue_to_alignment = {}
    residue_idx = 0
    col2res = {}
    col2score = {}

    for aln_idx, char in enumerate(aln):
        if char != '-':
            residue_idx += 1
            if residue_idx in important_residues:
                residue_to_alignment[residue_idx] = aln_idx
                col2res[aln_idx] = f"{aln[aln_idx]}{residue_idx}"
                col2score[aln_idx] = residue_to_score[residue_idx]
                
            # early stop if we've found everything
            if len(residue_to_alignment) == len(important_residues):
                break

    # print the corresponding residues in the original sequence
    if p:
        for residue, aln_idx in residue_to_alignment.items():
            print(f"Residue {residue} corresponds to alignment index {aln_idx}: {aln[aln_idx]}", flush=True)

    # store in dataframe        
    acc = [result.id for result in alignment]
    seqs = [list(str(result.seq)) for result in alignment]
    hits = [result.description.split(' ')[-1] for result in alignment]

    pssm = pd.DataFrame(seqs, index = acc)

    pssm = pssm.iloc[:, list(residue_to_alignment.values())]
    pssm['hit'] = hits
    pssm['contig'] = pssm['hit'].str.split('_').str[:-1].str.join('_')


    # check if cols contain correct residue for function
    def check_res(row):
        score = 0
        for col in residue_to_alignment.values():
            if row[col] == aln[col]:
                score += 1 + col2score[col]
            elif row[col] != '-': # higher weight for correct C
                    score += 1
        if score >= passing_score:
            return score
        else:
            return np.nan

    pssm['score'] = pssm.apply(check_res, axis = 1)
    pssm.dropna(subset = ['score'], inplace = True)
    pssm.rename(columns=col2res, inplace=True)

    return pssm

def check_gene_wrapper(gene_list: dict, config: dict, results_dir: str = "../results/tmp") -> dict:
    print("Checking for conserved residues...", flush=True)

    all_genes_checked = {}

    for gene in gene_list['Gene'].unique():
        print(f'checking {gene}', flush=True)

        # load json parameters
        ref_seq = config[gene]['ref_seq']
        important_residues = config[gene]['residues']
        residue_scores = config[gene].get('residue_scores', config[gene].get('reside_scores'))
        if residue_scores is None:
            residue_scores = [1] * len(important_residues)
        passing_score = config[gene]['passing_score']

        # initialize dataframe from first split
        first_file = os.path.join(results_dir, 'fasta_splits', f"{gene}_split.00001.aln")
        if not os.path.exists(first_file):
            all_genes_checked[gene] = pd.DataFrame()
            continue
        checked = check_gene(first_file, ref_seq, important_residues, residue_scores, passing_score, p=True)

        for file in glob.glob(os.path.join(results_dir, 'fasta_splits', f"{gene}_split.*.aln")):
            if f'{gene}_split.00001' in file:
                continue
            new = check_gene(file, ref_seq, important_residues, residue_scores, passing_score)
            checked = pd.concat([checked, new])

        checked.drop_duplicates(subset = ['hit'], inplace = True)
        checked.set_index(['hit'], append= True, inplace = True)
        checked.to_csv(os.path.join(results_dir, f'{gene}_residues.csv'))

        # only keep gene where residue matching passed
        tmp = gene_list[(gene_list['Gene'] == gene) & (gene_list['Hit'].isin(checked.index.get_level_values(1)))].copy()
        tmp['residue_match'] = gene
        tmp.to_csv(os.path.join(results_dir, f'{gene}_rescheck_nofilt.csv'))

        # potentially add code for (1) keeping best hit per genome and (2) removing hits that are in other genes (e.g. remove nifD,E hits from nifK)

        all_genes_checked[gene] = tmp

        print(f"{len(checked)} {gene} seqs \n", flush=True)

    # export
    os.makedirs(results_dir, exist_ok=True)
    all_genes_checked = pd.concat(all_genes_checked.values())
    all_genes_checked.to_csv(os.path.join(results_dir, 'nif_rescheck_nofilt.csv')) # index = GenomeID

    return all_genes_checked

# EXPORT
def export_results(
    all_genes_checked,
    min_genes: int,
    results_dir: str = "../results/final",
    proteins_dir: str = "../protein_faa_reps_232",
    include_metadata: bool = True,
    external: bool = False
) -> None:
    os.makedirs(results_dir, exist_ok=True)
    hits = all_genes_checked.copy()
    hits.set_index(['Gene', 'Hit'], inplace=True, append=True)
    hits['cluster'] = 'nif'

    # make sure nif clusters have at least nifHDK (unless it is a nifB-only operon)
    # doesn't check anf or vnf clusters
    def gene_check(genes):
        def has_gene_letter(letter):
            return any(letter in str(gene) for gene in genes)
        if len(genes) >= 3:
            if has_gene_letter('G'):
                return True
            elif has_gene_letter('H'):
                if has_gene_letter('D'):
                    if has_gene_letter('K'):
                        return True
        elif len(genes) == 1:
            if genes[0] == 'nifB':
                return True

        return False

    def relabel_operon_gene_prefix(idx, new_prefix):
        current = hits.loc[idx, 'residue_match'].astype(str)
        replace_mask = current.str.startswith(('nif', 'vnf', 'anf'))
        update_idx = current[replace_mask].index
        hits.loc[update_idx, 'residue_match'] = new_prefix + current.loc[update_idx].str[3:]
        hits.loc[update_idx, 'cluster'] = new_prefix

    for _, strand in hits.groupby(['GenomeID', 'contig', 'operon']):
        # Keep nifB unchanged while relabeling all other genes in the operon.
        non_nifb_mask = strand.index.get_level_values(1) != 'nifB'
        target_idx = strand.index[non_nifb_mask]

        # update to vnf if vnfG present, then update to anf if anfO present
        if strand['residue_match'].isin(['anfG', 'vnfG']).any():
            relabel_operon_gene_prefix(target_idx, 'vnf')

            # then update to anf if anfO present
            if strand['residue_match'].eq('anfO').any():
                relabel_operon_gene_prefix(target_idx, 'anf')

    hits.reset_index(inplace=True)
    # assign final annotation + group domains into proteins
    hits['protein'] = hits['Hit'].str.rsplit('_', n=1).str[0] # remove domain suffix to get protein name
    genomes_to_keep = pd.DataFrame(columns=hits.columns).set_index('protein')
    # assign annotation from remaining filtered seqs (now apply best hits approach)
    for idx, cluster in hits.groupby(['GenomeID', 'contig', 'operon', 'cluster']):
        # keep best hit per domain (of filtered domains)
        cluster = cluster.loc[cluster.groupby('Hit')['E-value'].idxmin()]
        # merge domains into protein, aggegrate rows
            # expecting nifNB, nifEN, nifDG
        if cluster['Orientation'].sum() < 0: # if all genes on negative strand, reverse order of genes
            cluster.sort_values(by=['Hit','Domain Start'], ascending=False, inplace=True)
        else:
            cluster.sort_values(by=['Hit', 'Domain Start'], ascending=True, inplace=True)
        # if 'LKPX01000013.1_53' in cluster['protein'].to_list():
        #     print(cluster)
        cluster = cluster.groupby('protein', group_keys=False).agg({
            'GenomeID': 'first',
            'contig': 'first',
            'operon': 'first',
            'cluster': 'first',
            'Gene': lambda x: '_'.join(sorted(x)),
            'Hit': lambda x: ','.join(sorted(x)),
            'E-value': 'min',
            'Bit Score': 'sum',
            'Bias': 'sum',
            'Domain Start': 'min',
            'Domain End': 'max',
            'Location': 'first',
            'Orientation': 'first',
            'Flag_Bias': 'sum',
            'pos_num': 'first',
            'top_hit': lambda x: ','.join(sorted(x)),
            'residue_match':  lambda x: ''.join(
                gene if i == 0 else gene[-1]
                for i, gene in enumerate(sorted(x))
            ),
            'GTDB': 'first'
        })

        # check if all nif genes are present (skip anf, vnf)
        if gene_check(cluster.Gene.to_list()):
            if genomes_to_keep.empty:
                genomes_to_keep = cluster.copy()
            else:
                genomes_to_keep = pd.concat([genomes_to_keep, cluster])

    #clean up cols
    genomes_to_keep['gene'] = genomes_to_keep['residue_match']
    genomes_to_keep.drop_duplicates(inplace = True)
    genomes_to_keep = genomes_to_keep.reset_index()

    # remove operons with fewer than 3 unique nif genes (except for operons that contain nifB)
    genomes_to_keep = filter_groups_by_unique_counts(
    genomes_to_keep,
    group_cols=["GenomeID", "contig", "operon"],
    requirements={
    "Gene": min_genes,
    "Hit": min_genes
    },
    exclude='nifB'
    )
    genomes_to_keep.set_index('GenomeID', inplace=True)

    # export csv with each gene as individual row
    print(external, flush=True)
    if external: # clean up export for external DiazoDB users
        external_export = genomes_to_keep.copy()
        external_export['Sequence Length'] = external_export['Domain End'] - external_export['Domain Start']
        external_export = external_export[['contig', 'protein', 'cluster', 'gene', 'Location', 'Orientation', 'Sequence Length']]
        external_export.to_csv(os.path.join(results_dir, 'nif_final.csv'))

    else:
        genomes_to_keep.to_csv(os.path.join(results_dir, 'nif_final.csv'))


    # get fasta sequences for each gene & export to fasta
    os.makedirs(os.path.join(results_dir, "fastas"), exist_ok=True)
    for gene in genomes_to_keep['gene'].unique().tolist():
        records = []

        for genome, hit in genomes_to_keep[genomes_to_keep.gene == gene].iterrows():
            record = get_seq(
                genome,
                hit.protein,
                id=hit.protein,
                description=f"{genome} {hit.gene}",
                dir=proteins_dir,
            )
            records.append(record)
                    
        # Write the records to a FASTA file
        with open(os.path.join(results_dir, "fastas/final_" + gene + ".fasta"), "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

        print(f"{len(records)} final {gene} hits", flush=True)

    # export csv grouped by genome, contig, operon, and cluster
    operons = genomes_to_keep.copy().set_index(['contig', 'operon', 'cluster'], append=True)
    operons.sort_values(by=['GenomeID', 'contig', 'operon', 'cluster', 'Location'], inplace=True)
    # Split Location by '-' and extract min/max of first and second parts
    location_parts = operons['Location'].str.split('-', expand=True)
    operons['_location_start'] = location_parts[0].astype(int)
    operons['_location_end'] = location_parts[1].astype(int)
    operons['Location_start'] = operons.groupby(['GenomeID', 'contig', 'operon', 'cluster'])['_location_start'].transform('min')
    operons['Location_end'] = operons.groupby(['GenomeID', 'contig', 'operon', 'cluster'])['_location_end'].transform('max')
    operons = operons.drop(columns=['_location_start', '_location_end'])

    operons = operons.groupby(['GenomeID', 'contig', 'operon', 'cluster']).agg(**{
        'Gene set': ('Gene', lambda x: ', '.join(x.astype(str).str.replace(r'^(nif|anf|vnf)', '', regex=True))),
        'Location_start': ('Location_start', 'first'),
        'Location_end': ('Location_end', 'last'),
        'GTDB': ('GTDB', 'first'),
        'Orientation': ('Orientation', 'first'),
        'pos_num': ('pos_num', lambda x: x.astype(int).tolist())
    })

    # def get_hit(rec):
    #     return rec.description.split(' ')[-1]

    # def get_cluster(file):
    #     clusters = {}
    #     clusters_fasta = list(SeqIO.parse(file, 'fasta'))

    #     index = 0
    #     while index < len(clusters_fasta)-1:
    #         rec = clusters_fasta[index]
    #         if rec.seq == '': # cluster header
    #             cluster = get_hit(clusters_fasta[index+1]) # hit
    #             clause = True
    #             members = []
    #             i = 1
    #             while clause & (index+i < len(clusters_fasta)): # add members (until next cluster header)
    #                 if clusters_fasta[index+i].seq != '':
    #                     members.append(get_hit(clusters_fasta[index+i]))
    #                     i+=1
    #                 elif clusters_fasta[index+i].seq == '':
    #                     clause = False
    #             index += i

    #             clusters[cluster] = members

    #     return clusters

    # operons['Group'] = ''
    # #     #genes = {'H': 'H', 'D': 'D_noOut'}
    # genes = {'H': 'H'}

    # for gene, file in genes.items():
    #     # get clustered datapoints
    #     clusters = get_cluster(f'../trees/nif{file}/clustered_nif{file}_all_seqs.fasta') 

    #     # assign group 
    #     for group in ['1', '2', '3', '4a', '4c', '3anfvnf']:
    #         lines = []
    #         hits = []
    #         with open(f'nif_groups/nif{gene}_group{group}.txt','r') as f:
    #             lines = f.read().splitlines()
    #             for line in lines:
    #                 hit = '_'.join(line.split('|')[-1].split(' '))[:-1]
    #                 hits.append(hit) # reformat "hits" to match nif index
    #                 hits.extend(clusters[hit]) # add clustered hits to list of hits to update
    #         for hit in hits:
    #             operons.loc[operons.index == hit, 'Group'] = f'Group {group}'

#     # group by genome, contig and save
#     nif = nif[['GenomeID', 'contig', 'Gene set', 'Position', 'GTDB', 'Location' ,'Orientation', 'Group']]
#     nif = nif.groupby(['GenomeID','contig']).first()
#     nif.reset_index(level=['GenomeID','contig'], inplace = True)

    if include_metadata:
        # Database-build mode enriches representative genomes with taxonomy.
        gtdb_metadata = pd.read_csv(
            "GTDB_metadata.gz",
            sep="\t",
            usecols=[
                "accession",
                "gtdb_taxonomy",
                "ncbi_taxonomy",
                "ncbi_taxonomy_unfiltered",
                "ncbi_country",
                "ncbi_isolation_source",
            ],
        )

        operons = operons.join(
            gtdb_metadata.set_index("accession"),
            on="GenomeID",
            how="left",
        ).drop(columns=["GTDB"])

        operons["Organism"] = (
            operons["gtdb_taxonomy"].str.split(";").str[-1].str.split("__").str[-1]
        )
        operons["Regulon"] = ""
        operons["Operon"] = ""
        operons["Group"] = ""
        operons.rename(
            columns={
                "Gene set": "Nitrogenase Set",
                "Group": "Group No",
                "gtdb_taxonomy": "GTDB Taxonomy",
                "ncbi_taxonomy": "NCBI Taxonomy",
                "ncbi_isolation_source": "Isolation Source",
            },
            inplace=True,
        )
    else:
        # Uploaded genomes do not need database taxonomy enrichment.
        operons.drop(columns=["GTDB"], errors="ignore", inplace=True)
        operons.rename(columns={"Gene set": "Nitrogenase Set"}, inplace=True)

    if external: # clean up export for external DiazoDB users
        # append contig to the front of each value in pos_num 
        operons['proteins'] = operons.apply(lambda row: [f"{row['contig']}_{pos}" for pos in row['pos_num']], axis=1)
        operons = operons[['contig', 'proteins', 'cluster', 'Nitrogenase Set', 'Location_start', 'Location_end', 'Orientation']]

    operons.to_csv(os.path.join(results_dir, 'nif_clusters.csv'))

    return genomes_to_keep

def main() -> None:
    args = parse_args()

    # load config and nif table at runtime
    config_file = args.config_file
    config = json.load(open(config_file, 'r'))

    # load hits dataframe
    hits = load_nif(update_index=['GenomeID'], hits_file=args.hits_file)
    
    if hits.empty:
        print("No qualifying HMM hits; skipping conserved residue check.", flush=True)
        return

    # store per-gene df in a dict
    gene_list = {gene: hits[hits.Gene == gene].copy() for gene in config.keys()}
    print(f"checking conserved residues for: {list(gene_list)}", flush=True)

    # reload + align fasta files if specified
    if args.reload_fasta:
        gene_list = write_fasta(gene_list, config, proteins_dir=args.proteins_dir, results_dir=args.results_dir)
        args.align_fasta = True  # Automatically align if reloading fasta
    else: 
        gene_list = pd.read_csv(os.path.join(args.results_dir, 'hits_seq.csv'))
    if args.align_fasta:
        aln_fasta(gene_list, config, results_dir=args.results_dir)

    # check for conserved residues and save results
    checked = check_gene_wrapper(gene_list, config, results_dir=args.results_dir)
    export_results(
        checked,
        args.min_genes,
        results_dir=args.final_dir,
        proteins_dir=args.proteins_dir,
        include_metadata=not args.skip_metadata,
        external=args.external,
    )

    # check results against known nitrogenase
    if args.test:
        test(checked=pd.read_csv('../results/final/nif_final.csv'))

if __name__ == "__main__":
    main()
