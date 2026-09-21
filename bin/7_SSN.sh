#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=2:15:00   # walltime # about 3hrs for ~300 seqs, 20+ hrs for 7000+
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 10GB   # memory per node
#SBATCH --job-name=SSN   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# load env, software
eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/make_trees

module load mafft/7.505-gcc-13.2.0-nklkvtc

DIR="../diazoDB-comparison/SSN"
CLUSTER="${DIR}/nifD_anfD_vnfD_12192023.faa"

# cluster NFixDB @ 0.85% ID
mmseqs createdb "$CLUSTER" tmp/seqDB
mmseqs cluster tmp/seqDB tmp/clustered tmp --min-seq-id 0.85 -c 0.8 --cov-mode 0
mmseqs createtsv tmp/seqDB tmp/seqDB tmp/clustered "${CLUSTER}.tsv"
mmseqs result2repseq tmp/seqDB tmp/clustered tmp/clustered_reps
mmseqs result2flat tmp/seqDB tmp/seqDB tmp/clustered_reps "$CLUSTER" --use-fasta-header

# add DB prefix
sed 's/^>/&NFixDB_/' "${DIR}/nifD_anfD_vnfD_12192023_clustered.fasta" > "${DIR}/nifD_anfD_vnfD_12192023_clustered.fasta"
sed 's/^>/&NSDB_/' "${DIR}/converted-nifD-extant.fasta" > "${DIR}/converted-nifD-extant.fasta"
sed 's/^>/&Nif-Finder_/' "${DIR}/true-nifD.faa" > "${DIR}/true-nifD.faa"
sed 's/^>/&DiazoDB_/' "${DIR}/nifD_anfD_vnfD_clustered.fasta" > "${DIR}/nifD_anfD_vnfD_clustered.fasta"

# concat all DB seqs
cat "${DIR}/converted-nifD-extant.fasta" "${DIR}/nifD_anfD_vnfD_12192023_clustered.fasta"  "${DIR}/nifD_anfD_vnfD_clustered.fasta"  "${DIR}/true-nifD.faa" > "${DIR}/nifD-comparison.fasta"

# blastn (all-against-all)
mmseqs createdb "${DIR}/nifD-comparison.fasta" "${DIR}/queryDB"
mmseqs createdb "${DIR}/nifD-comparison.fasta" "${DIR}/targetDB"
mmseqs createindex "${DIR}/targetDB" "${DIR}/tmp"
mmseqs search "${DIR}/queryDB" "${DIR}/targetDB" "${DIR}/resultDB" tmp
mmseqs convertalis "${DIR}/queryDB" "${DIR}/targetDB" "${DIR}/resultDB" "${DIR}/resultDB.m8"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
