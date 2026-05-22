"""Parse HMM output and store significant hits.

Inputs:
    --hits    Glob or path to hmmsearch.out results from the HMM search step.
    --outdir  Output directory for hits.feather and hits.csv.

Flags:
    --reload-fasta  Accepted for pipeline compatibility; not used in this step.

Example:
    python bin/Parse_hmm_results.py \
        --hits "../results/hmmsearch_results/archaea/*.out" \
        --outdir ../results/archaea

Hits must meet the following criteria:
    - positive bit score
    - full sequence evalue < 0.01
    - best domain evalue < 0.01 (otherwise flagged for review)
"""


import pandas as pd
import re
import glob
import argparse
from pathlib import Path


from Bio import SearchIO

# get outpur folder
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse hmmsearch outputs and write hits.feather and hits.csv."
        )
    )
    parser.add_argument(
        "--hits",
        required=True,
        help="Path to hmmsearch.out from the HMM search step.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for hits.feather and hits.csv.",
    )
    parser.add_argument(
        "--reload-fasta",
        action="store_true",
        help=(
            "Reload FASTA sequences for downstream steps if needed. Lengthy process so skip if not neccessary."
            "This flag is accepted for pipeline compatibility but is not used here."
        ),
    )
    return parser.parse_args()

# Get taxonomy from GTDB, https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/
    # NCBI taxonomy from metadata files

def parse_hits(hits_path: Path, outdir: Path) -> None:

    GTDB_taxonomy = pd.read_csv('GTDB_taxonomy.tsv', header = None, sep = '\t', names=['GenomeID', 'GTDB'])

    # Create empty lists
    result_target = []
    query_id = []
    hit_id = []
    evalue = []
    best_domain_evalue = []
    bitscore = []
    bias = []
    location = []
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
        str = r'# ([0-9]+) # ([0-9]+)'
        location.append(re.match(str, item.description).group(1) + "-" + re.match(str, item.description).group(2))
        # grab full alignment length (need to sum all domains)
        alen = 0
        for domain in item.hsps:
            alen += domain.aln_span
        alens.append(alen)
        slength.append(int(re.match(str, item.description).group(2))-int(re.match(str, item.description).group(1)))

    # Parse through files in output directory

    for file in glob.glob(str(hits_path)):

        # RegEx for the GenomeID (double checking that file is really a genome)
        try:
            str = r'([\w]+_[\w]+_[\d]+\.[\d])'
            genomeID = re.search(str, file).group()
        except:
            continue

        # Parse file using SearchIO/HmmerIO
        for result in SearchIO.parse(file, 'hmmer3-text'):
            for item in result.hits:

                # grab gene name
                str = r'([a-zA-Z]+)' # ex. nifHDK, pchlide
                gene = re.findall(str, result.id)[0]

                # Check for positive bitscore and append the data to the corresponding lists
                if item.bitscore > 0 and item.evalue < 0.01:
                    # append hits
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
                    elif item.bias == 0:
                        flag2.append(0)
                    else:
                        flag2.append(1)
                    
    # create and store dataframe
    hits = pd.DataFrame({'GenomeID': result_target, 'Gene': query_id, 'Hit': hit_id, 
                        'E-value': evalue, 'Best Domain E-value': best_domain_evalue, 'Bit Score': bitscore, 'Bias': bias,
                        'Location': location, 'Alignment Length': alens, 'Sequence Length': slength, 
                        'Flag_Eval': flag1, 'Flag_Bias': flag2})

    # add taxonomy info
    hits = pd.merge(hits, GTDB_taxonomy, on = "GenomeID", how = "left")

    # save as feather file
    outdir.mkdir(parents=True, exist_ok=True)
    hits.to_feather(outdir / "hits.feather")
    hits.to_csv(outdir / "hits.csv", index=False)

def main() -> None:
    args = parse_args()
    parse_hits(Path(args.hits), Path(args.outdir))


if __name__ == "__main__":
    main()