#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=3:15:00   # walltime # about 3hrs for ~300 seqs, 20+ hrs for 7000+
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per node
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

# add outgroup
#cat ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_09052025.fasta ../trees/arsA_GTDB_WW_GET3_09052025/GET3-sampled.fasta > ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.fasta

# cluster
#mmseqs easy-cluster ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.fasta ../trees/arsA_GTDB_WW_GET3_09052025/clustered_arsA_GTDB_WW_GET3 tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0

## count clusters
#num=$(grep ">" ../trees/arsA_GTDB_WW_GET3_09052025/clustered_arsA_GTDB_WW_GET3_rep_seq.fasta | wc -l)
#echo "$num clusters for 0.9"

# align ars sequences
#mafft --auto --thread 4 ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.fasta > ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.aln

# remove gappy alignments
#trimal -in ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.aln -out ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.trim -sgc -gappyout -keepheader

echo "tree building"
# build maximum likelihood tree
iqtree -s ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000 #use this to find best model and threads
#iqtree -s ../trees/arsA_GTDB_WW_GET3_09052025/arsA_GTDB_WW_GET3_09052025.trim -safe -m MFP -msub nuclear -T 8 -B 1000 -alrt 1000 #run w correct threads and model

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
