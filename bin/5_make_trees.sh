#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=20:05:00   # walltime # about 3hrs for ~300 seqs, 20+ hrs for 7000+
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 30GB   # memory per node
#SBATCH -J tree.e%j   # job name

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=zshivji@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## /SBATCH -p general # partition (queue)
## /SBATCH -o slurm.%N.%j.out # STDOUT
## /SBATCH -e slurm.%N.%j.err # STDERR

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
conda activate /home/zshivji/miniconda3/envs/make_trees/

module load mafft/7.505-gcc-13.2.0-nklkvtc

echo "preprocessing"
# cluster
mmseqs easy-cluster ../results/checked_nifK04292025.fasta ../trees/nifK_noOut_04292025/clustered_nifK_noOut tmp --min-seq-id 0.74 -c 0.8 --cov-mode 0

# count clusters
num=$(grep ">" ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.fasta | wc -l)
echo "$num clusters for 0.77"

# add outgroup
#cat ../trees/BchB.fasta >> ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.fasta

# align nif sequences
mafft --auto --thread 4 ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.fasta > ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.aln

# remove gappy alignments
trimal -in ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.aln -out ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -sgc -gappyout -keepheader

echo "tree building"
# build maximum likelihood tree
#iqtree -s ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000 #use this to find best model and t>
iqtree -s ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -safe -m LG+R10 -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000 #use this to find best model and threads

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
