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

from helper import load_nif, filter_groups_by_unique_counts, get_seq
from test import test


def default_proteins_dir() -> str:
    """Resolve the GTDB protein reps directory with sensible fallbacks."""
    if os.getenv("DIAZODB_PROTEIN_REPS_DIR"):
        return os.environ["DIAZODB_PROTEIN_REPS_DIR"]
    if os.path.isdir("../protein_faa_reps_latest"):
        return "../protein_faa_reps_latest"
    return "../protein_faa_reps_232"

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
            "Path to GTDB representative protein directories. "
            "Defaults to DIAZODB_PROTEIN_REPS_DIR, then ../protein_faa_reps_latest, "
            "then ../protein_faa_reps_232."
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

    return parser.parse_args()

# Note: runtime data (nif table, config, per-gene frames) are loaded in `main`.
# Helper functions below accept explicit inputs and return results.

def write_fasta(gene_list: dict, config: dict, proteins_dir: str = "../protein_faa_reps_232", results_dir: str = "../results/tmp"):
    """Get fasta sequences for each gene & export to fasta.

    Returns updated gene_list (with 'Seq' column filled where found).
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
            seq = get_seq(genome, hit.Hit[:-2], start=hit['Domain Start'], end=hit['Domain End'])
            records.append(seq)

        print(len(records), "records found for " + name, flush=True)
        with open(os.path.join(results_dir, name + ".fasta"), "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

    gene_list = pd.concat(gene_list.values())
    gene_list.to_csv(os.path.join(results_dir, 'hits_seq.csv'))

    return gene_list

def aln_fasta(gene_list: dict, results_dir: str = "../results/tmp"):
    # align fasta files
    splits_dir = os.path.join(results_dir, 'fasta_splits')
    os.makedirs(splits_dir, exist_ok=True)
    for name, gene in gene_list.items():
        print("aligning "+ name, flush=True)
        num = int(gene.shape[0]/200)+1 # how many splits
        os.makedirs("../results/tmp/fasta_splits", exist_ok=True)
        os.system(f"seqtk split -n {num} ../results/tmp/fasta_splits/{name}_split ../results/tmp/{name}.fasta") # split fasta file
        for i in range(num):
            os.system(f"seqtk subseq ../results/tmp/{name}.fasta ../results/ref_seq.ids >> ../results/tmp/fasta_splits/{name}_split.{str(i+1).zfill(5)}.fa") # add reference sequences
            os.system(f"mafft --auto --quiet --thread 4 ../results/tmp/fasta_splits/{name}_split.{str(i+1).zfill(5)}.fa > ../results/tmp/fasta_splits/{name}_split.{str(i+1).zfill(5)}.aln")


def check_gene(file, ref_seq, important_residues, residue_scores, passing_score, p=False):
    if len(important_residues) != len(residue_scores):
        raise ValueError(
            f"{file[-20:-16]}: residues and residue_scores must have the same length "
            f"({len(important_residues)} != {len(residue_scores)})"
        )

    residue_to_score = dict(zip(important_residues, residue_scores))
    alignment = AlignIO.read(f"{file}", "fasta")

    for record in alignment:
        if ref_seq in record.description:
            aln = record.seq
            break

    # map important cols in aln to ref_seq
    residue_to_alignment = {}
    residue_idx = 0  # index in original (ungapped) sequence
    col2res = {}
    col2score = {}

    for aln_idx, char in enumerate(aln):
        if char != '-':
            residue_idx += 1
            if residue_idx in important_residues:
                residue_to_alignment[residue_idx] = aln_idx
                col2res[aln_idx] = aln[aln_idx]
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
        ref_seq = config[gene]['ref_gene']
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
    all_genes_checked.to_csv(os.path.join(results_dir, 'nif_rescheck_nofilt.csv'), index=False)

    return all_genes_checked

# EXPORT
def export_results(all_genes_checked, min_genes: int, results_dir: str = "../results/final") -> None:
    hits = all_genes_checked.copy()
    hits.set_index(['GenomeID', 'Gene', 'Hit'], inplace=True)
    hits['cluster'] = 'nif'

    # make sure nif clusters have at least nifHDK (unless it is a nifB-only operon)
    # doesn't check anf or vnf clusters
    def gene_check(genes):
        if len(genes) >= 3:
            if genes.__contains__('vnfG'):
                return True
            elif genes.__contains__('anfO'):
                return True
            elif genes.__contains__('nifH'):
                if genes.__contains__('nifD'):
                    if genes.__contains__('nifK'):
                        return True
        else:
            if genes.__contains__('nifB'):
                return True

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

    # remove operons with fewer than 3 unique nif genes (except for operons that contain nifB)
    hits.reset_index(inplace=True)
    hits = filter_groups_by_unique_counts(
        hits,
        group_cols=["GenomeID", "contig", "operon"],
        requirements={
        "Gene": min_genes,
        "Hit": min_genes
        },
        exclude='nifB'
    )

############# Need to get combined Gene and Seq Location? fixed ##############
    # insted of Domain start/end min/max --> first/last now that it's sorted coretly (based on strand +/-1)

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
    genomes_to_keep['Gene'] = genomes_to_keep['residue_match']
    #hits = hits[['Gene', 'E-value', 'Bit Score', 'Location', 'Orientation', 'Alignment Length', 'Sequence Length', 'GTDB']]
    genomes_to_keep.drop_duplicates(inplace = True)
    genomes_to_keep = genomes_to_keep.reset_index().set_index('GenomeID')

    # export csv with each gene as individual row
    genomes_to_keep.to_csv(os.path.join(results_dir, 'nif_final.csv'))

    # get fasta sequences for each gene & export to fasta
    os.makedirs(os.path.join(results_dir, "fastas"), exist_ok=True)
    for gene in genomes_to_keep['Gene'].unique().tolist():
        records = []

        for genome, hit in genomes_to_keep[genomes_to_keep.Gene == gene].iterrows():
            record = get_seq(genome, hit.protein, id=hit.protein, description=f"{genome} {hit.Gene}")
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
        'Gene set': ('Gene', lambda x: ', '.join(x.astype(str))),
        'Location_start': ('Location_start', 'first'),
        'Location_end': ('Location_end', 'last'),
        'GTDB': ('GTDB', 'first'),
        'Orientation': ('Orientation', 'first')
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

    #merge dataframes to include accession, metadata, and sequences (for filtering)
    GTDB_metadata = pd.read_csv('GTDB_metadata.gz', sep = '\t', usecols=['accession', 'gtdb_taxonomy', 'ncbi_taxonomy', 
                                                    'ncbi_taxonomy_unfiltered', 'ncbi_country', 'ncbi_isolation_source'])
    
    operons = operons.join(
        GTDB_metadata.set_index('accession'),
        on='GenomeID',
        how='left'
    ).drop(columns=['GTDB'])

    operons['Organism'] = operons['gtdb_taxonomy'].str.split(';').str[-1].str.split('__').str[-1]
    operons['Regulon'] = ''
    operons['Operon'] = ''
    operons['Group'] = ''
#    nif['PredGrowthTemp'] = ''
    operons.rename(columns = {'GenomeID': 'Genome', 'contig': 'Contig', 'Gene set': 'Nitrogenase Set', 
                        'Group':'Group No', 'gtdb_taxonomy': 'GTDB Taxonomy', 'cluster': 'Cluster', 
                        'ncbi_taxonomy': 'NCBI Taxonomy', 'ncbi_isolation_source':'Isolation Source'}, inplace = True)

    operons.to_csv(f'../results/final/nif_clusters.csv')

    return genomes_to_keep

#----------------------------------------------------------------

    # for gene_name, checked_df in all_genes_checked.items():
    #     for genome_idx, cols in checked_df.iterrows():
    #         try:
    #             nif.loc[(genome_idx[0], genome_idx[1]), 'residue_match'] = gene_name
    #         except Exception:
    #             continue

    # # filter to get hits that passed residue matching and are within length bounds
    # filtered = {}
    # for gene_name, params in config.items():
    #     maxl = params.get('max-length')
    #     minl = params.get('min-length')
    #     mask = (nif.residue_match == gene_name) & (nif.Gene == gene_name)
    #     if minl is not None:
    #         mask &= (nif['Alignment Length'] >= minl)
    #     if maxl is not None:
    #         mask &= (nif['Alignment Length'] < maxl)
    #     filtered[gene_name] = nif[mask]


def main() -> None:
    args = parse_args()

    # load config and nif table at runtime
    config_file = args.config_file
    config = json.load(open(config_file, 'r'))

    # load hits dataframe
    hits = load_nif(update_index=['GenomeID'])
    if 'Seq' not in hits.columns:
        hits['Seq'] = ''

    # store per-gene df in a dict
    gene_list = {gene: hits[hits.Gene == gene].copy() for gene in config.keys()}
    print(f"checking conserved residues for: {list(gene_list)}", flush=True)

    # reload + align fasta files if specified
    if args.reload_fasta:
        gene_list = write_fasta(gene_list, config, proteins_dir=args.proteins_dir)
        args.align_fasta = True  # Automatically align if reloading fasta
    else: 
        gene_list = pd.read_csv(os.path.join("../results/tmp", 'hits_seq.csv'))
    if args.align_fasta:
        aln_fasta(gene_list)

    # check for conserved residues and save results
    checked = check_gene_wrapper(gene_list, config)
    export_results(checked, args.min_genes)

    # check results against known nitrogenase
    if args.test:
        test(checked)

if __name__ == "__main__":
    main()



#### END OF CONSERVED RESIDUE CHECKING, BELOW IS VALIDATion/BACKUP CODE ####

# #backup check --> likely don't need this since we're reducing alignments to 200 seq per file
# #get seq that failed checks
# print('getting failed sequences', flush=True)

# seqs = []

# for gene in 'DKEN':

#     # for each gene, get all hmm hit acc
#     result = list(SeqIO.parse(f"../results/fasta_splits/nif{gene}.fasta", "fasta"))
#     hit = [record.description.split(" ")[-1] for record in result]

#     # get seq that failed check
#     for record, acc in zip(result, hit):
#         if acc not in eval(f"list(nif{gene}_checked.index.get_level_values(1).unique())"):
#             seq = SeqRecord(Seq(record.seq), id=record.id, description=acc)
#             seqs.append(record)

# print(str(len(seqs)) + " seqs failed nifDKEN checks", flush=True) # counting how many failed pre/post 2026 edits

# # Write the records to a FASTA file
# with open(f"../results/fasta_splits/nifDKEN.fasta", "w") as output_handle:
#     SeqIO.write(seqs, output_handle, "fasta")

# # add reference sequences
# for gene in 'DKEN':
#     os.system(f"seqtk subseq ../results/fasta_splits/nif{gene}.fasta ../results/ref_seq.ids >> ../results/fasta_splits/nifDKEN.fasta")

# print('aligning failed sequences', flush=True)
# # aln all seqs
# num = int(len(seqs)/200) +1  # how many splits
# os.system(f"seqtk split -n {num} ../results/fasta_splits/nifDKEN_split ../results/fasta_splits/nifDKEN.fasta") # split fasta file
# for i in range(num):
#     print(i+1)
#     os.system(f"seqtk subseq ../results/fasta_splits/nifDKEN.fasta ../results/ref_seq.ids >> ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.fa") # add reference sequences
#     os.system(f"mafft --auto --quiet --thread 4 ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.fa > ../results/fasta_splits/nifDKEN_split.{str(i+1).zfill(5)}.aln")
    

# backup check --> likely don't need this since we're reducing alignments to 200 seq per file
# print('checking failed sequences', flush=True)

# nifD_backup = check_gene('nifDKEN_split.00001', ref_seq_nifD, important_residues_nifD, passing_score)
# nifK_backup = check_gene('nifDKEN_split.00001', ref_seq_nifK, important_residues_nifK, passing_score)
# nifE_backup = check_gene('nifDKEN_split.00001', ref_seq_nifE, important_residues_nifE, passing_score)
# nifN_backup = check_gene('nifDKEN_split.00001', ref_seq_nifN, important_residues_nifN, passing_score)

# for file in glob.glob(f'nifDKEN_split.00*.aln'):
#     if file == '../results/fasta_splits/nifDKEN_split.00001.aln':
#         continue
#     nifD_backup = check_gene(file[:-4], ref_seq_nifD, important_residues_nifD, passing_score)
#     nifD_backup = pd.concat([nifD_backup, new])

#     nifK_backup = check_gene(file[:-4], ref_seq_nifK, important_residues_nifK, passing_score)
#     nifK_backup = pd.concat([nifK_backup, new])

#     nifE_backup = check_gene(file[:-4], ref_seq_nifE, important_residues_nifE, passing_score)
#     nifE_backup = pd.concat([nifE_backup, new])

#     nifN_backup = check_gene(file[:-4], ref_seq_nifN, important_residues_nifN, passing_score)
#     nifN_backup = pd.concat([nifN_backup, new])

# # set index as genome, hit
# nifD_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifD_backup = nifD_backup[~nifD_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per contig
# nifD_backup['gene_cluster'] = 0
# nifD_backup = nifD_backup.loc[nifD_backup.groupby(['contig'])['score'].idxmax()]
# nifD_backup.drop_duplicates(inplace = True)

# print(str(len(nifD_backup.index.unique())) + " nifD seqs", flush=True)
# nifD_checked = pd.concat([nifD_checked, nifD_backup])

# # set index as genome, hit
# nifK_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifK_backup = nifK_backup[~nifK_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per contig
# nifK_backup['gene_cluster'] = 0
# nifK_backup = nifK_backup.loc[nifK_backup.groupby(['contig'])['score'].idxmax()]
# nifK_backup.drop_duplicates(inplace = True)

# print(str(len(nifK_backup.index.unique())) + " nifK seqs", flush=True)
# nifK_checked = pd.concat([nifK_checked, nifK_backup])

# # set index as genome, hit
# nifE_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifE_backup = nifE_backup[~nifE_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per gene cluster
# nifE_backup['gene_cluster'] = 0
# nifE_backup = nifE_backup.loc[nifE_backup.groupby(['contig'])['score'].idxmax()]
# nifE_backup.drop_duplicates(inplace = True)

# print(str(len(nifE_backup.index.unique())) + " nifE seqs", flush=True)
# nifE_checked = pd.concat([nifE_checked, nifE_backup])

# # set index as genome, hit
# nifN_backup.set_index(['hit'], append= True, inplace = True)

# # remove hits that are already in saved
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifD_checked.index.get_level_values(1).to_list())] # remove nifD
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifK_checked.index.get_level_values(1).to_list())] # remove nifK
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifE_checked.index.get_level_values(1).to_list())] # remove nifE
# nifN_backup = nifN_backup[~nifN_backup.index.get_level_values(1).isin(nifN_checked.index.get_level_values(1).to_list())] # remove nifN

# # for each genome, only keep the best hit per gene cluster
# nifN_backup['gene_cluster'] = 0
# nifN_backup = nifN_backup.loc[nifN_backup.groupby(['contig'])['score'].idxmax()]
# nifN_backup.drop_duplicates(inplace = True)

# print(str(len(nifN_backup.index.unique())) + " nifN seqs", flush=True)
# nifN_checked = pd.concat([nifN_checked, nifN_backup])

# # append updated annotation (based on conserved residue matching) to nif files
# nif['residue_match'] = ''
# nif['backup_match'] = ''

# nif.set_index(['Hit'], append = True, inplace = True)
# nif.sort_index(inplace = True)

# # update residue match column in nif df
# for gene in 'HDKEN':
#     for genome, cols in eval(f"nif{gene}_checked.iterrows()"):
#         nif.loc[(genome, cols.hit), 'residue_match'] = "nif" + gene

#  # add backup check       
# for gene in 'DKEN':
#     for genome, cols in eval(f"nif{gene}_backup.iterrows()"):
#         nif.loc[(genome[0], genome[1]), 'backup_match'] = "nif" + gene

# # filter to get hits that passed residue matching
# nifH = nif[(nif.residue_match == 'nifH') & (nif.Gene == 'nifH') & (nif['Alignment Length'] > 200)]
# # nifB = nif[(nif.Gene == 'nifB')] # not done
# # nifB['residue_match'] = 'nifB'

# # only index OG matches
# nifD = nif[((nif.residue_match == 'nifD') & (nif.Gene == 'nifD') & (nif.backup_match != 'nifD') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifD_backup = nif[(nif.backup_match == 'nifD') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifD_backup = nifD_backup.groupby(['GenomeID', 'Hit']).first()
# nifD = pd.concat([nifD, nifD_backup])

# # only index OG matches
# nifK = nif[((nif.residue_match == 'nifK') & (nif.Gene == 'nifK') & (nif.backup_match != 'nifK') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifK_backup = nif[(nif.backup_match == 'nifK') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifK_backup = nifK_backup.groupby(['GenomeID', 'Hit']).first()
# nifK = pd.concat([nifK, nifK_backup])

# # only index OG matches
# nifE = nif[((nif.residue_match == 'nifE') & (nif.Gene == 'nifE') & (nif.backup_match != 'nifE') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifE_backup = nif[(nif.backup_match == 'nifE') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifE_backup = nifE_backup.groupby(['GenomeID', 'Hit']).first()
# nifE = pd.concat([nifE, nifE_backup])

# # only index OG matches
# nifN = nif[((nif.residue_match == 'nifN') & (nif.Gene == 'nifN') & (nif.backup_match != 'nifN') & (nif['Alignment Length'] > 300))]
# # add backup check
# nifN_backup = nif[(nif.backup_match == 'nifN') & (nif['Alignment Length'] > 300)].sort_values(by = 'E-value')
# nifN_backup = nifN_backup.groupby(['GenomeID', 'Hit']).first()
# nifN = pd.concat([nifN, nifN_backup])

# nif = pd.concat([nifH, nifD, nifK, nifE, nifN])
# nif.sort_index(inplace = True)

# nif.to_csv(f'../results/fasta_splits/nif_rescheck_nofilt.csv')

# # Are any nifD,E,K being missed and printed nifN (for example?) should I align all DKEN first and then check for conserved residues?

# EXPORT
