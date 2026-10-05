#!/bin/bash

#SBATCH --time=01:04:00   # walltime #8hrs?
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH --job-name=checkDB   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

ROOT="/resnick/groups/enviromics/zahra/diazoDB-HPC"
CONFIG_FILE="nif-config.json"
DIR="${ROOT}/diazoDB-comparison"

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 <file-gene-table.csv>" >&2
    echo "CSV must have headers: file,gene" >&2
    exit 2
fi

TABLE="$1"
if [[ "${TABLE}" = /* ]]; then
    TABLE_PATH="${TABLE}"
elif [[ -e "${TABLE}" ]]; then
    TABLE_PATH="${TABLE}"
elif [[ -e "${ROOT}/${TABLE}" ]]; then
    TABLE_PATH="${ROOT}/${TABLE}"
elif [[ -e "${DIR}/${TABLE}" ]]; then
    TABLE_PATH="${DIR}/${TABLE}"
else
    echo "CSV table not found: ${TABLE}" >&2
    exit 1
fi

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : ${SLURM_JOBID:-NA} / ${SLURM_JOB_NAME:-checkDB}"
echo "======================================================"
echo ""

module load mafft/7.505-gcc-13.2.0-nklkvtc

tail -n +2 "${TABLE_PATH}" | while IFS=, read -r file gene; do
    file="${DIR}/${file}"
    file_dir="$(dirname "${file}")"

    echo "Processing ${file} as ${gene}"

    python - "${file}" "${file_dir}" "${CONFIG_FILE}" "${gene}" <<'PY'
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

from conserved_res import check_gene

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fasta = Path(sys.argv[1])
file_dir = Path(sys.argv[2])
config_file = Path(sys.argv[3])
gene = sys.argv[4]

config = json.load(open(config_file, 'r'))

splits_dir = file_dir / "fasta_splits"
os.makedirs(splits_dir, exist_ok=True)

num_records = sum(1 for _ in SeqIO.parse(fasta, "fasta"))
num_splits = int(num_records / 200) + 1
split_prefix = splits_dir / f"{fasta.stem}_split"

subprocess.run(
    ["seqtk", "split", "-n", str(num_splits), str(split_prefix), str(fasta)],
    check=True,
)

important_residues = config[gene]["residues"]
residue_scores = config[gene].get("residue_scores", config[gene].get("reside_scores"))
if residue_scores is None:
    residue_scores = [1] * len(important_residues)
passing_score = config[gene]["passing_score"]

checked = []
for i, split_path in enumerate(sorted(glob.glob(f"{split_prefix}.*.fa")), start=1):
    alignment_path = f"{split_path[:-3]}.aln"

    ref = SeqRecord(Seq(config[gene]["ref_seq"]), id="reference", description=gene)
    with open(split_path, "a") as f:
        SeqIO.write(ref, f, "fasta")

    with open(alignment_path, "wb") as alignment_file:
        subprocess.run(
            ["mafft", "--auto", "--quiet", "--thread", "4", split_path],
            stdout=alignment_file,
            check=True,
        )

    # append the results of each split, but discard the reference row (last row)
    checked.append(
        check_gene(
            alignment_path,
            important_residues,
            residue_scores,
            passing_score,
            p=(i == 1),
        )
    )

df = pd.concat(checked) if checked else pd.DataFrame()
df.drop_duplicates(inplace=True)
out_file = file_dir / f"{fasta.stem}_rescheck.csv"
df.to_csv(out_file)
print(f"Saved {df.shape[0]} rows to {out_file} \n", flush=True)
PY
done

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
