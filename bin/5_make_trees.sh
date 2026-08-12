#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=10:05:00   # walltime # about 3hrs for ~300 seqs, 20+ hrs for 7000+
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 30GB   # memory per node
#SBATCH --job-name=tree   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# load env, software
eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/make_trees

module load mafft/7.505-gcc-13.2.0-nklkvtc

echo "preprocessing"
# cluster, to keep full fasta header, run easy-cluster workflow separately
cat ../results/final/fastas/final_nifH.fasta ../results/final/fastas/final_vnfH.fasta ../results/final/fastas/final_anfH.fasta > ../trees/nifH/nifH_vnfH_anfH.fasta
mmseqs createdb ../trees/nifH/nifH_vnfH_anfH.fasta tmp/seqDB
mmseqs cluster tmp/seqDB tmp/clustered tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs result2repseq tmp/seqDB tmp/clustered tmp/clustered_reps
mmseqs result2flat tmp/seqDB tmp/seqDB tmp/clustered_reps ../trees/nifH/nifH_vnfH_anfH_clustered.fasta --use-fasta-header

# count clusters
num=$(grep ">" ../trees/nifH/nifH_vnfH_anfH_clustered.fasta | wc -l)
echo "$num clusters for 0.9"

# add outgroup
cat ../trees/BchL.fasta >> ../trees/nifH/nifH_vnfH_anfH_clustered.fasta


#{
#awk '/^>/{print ">nifD|" substr($0,2); next} {print}' ../HMM_seeds/clustered_nifD_rep_seq.fasta
#awk '/^>/{print ">nifK|" substr($0,2); next} {print}' ../HMM_seeds/clustered_nifK_rep_seq.fasta
#awk '/^>/{print ">nifE|" substr($0,2); next} {print}' ../HMM_seeds/clustered_nifE_rep_seq.fasta
#awk '/^>/{print ">nifN|" substr($0,2); next} {print}' ../HMM_seeds/nifN_merged_len300.fasta
#awk '/^>/{print ">nifB|" substr($0,2); next} {print}' ../HMM_seeds/nifB.fasta
#} > ../HMM_seeds/nifDKENB.fasta

# align nif sequences
#mafft --auto --thread 4 ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.fasta > ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.aln
#mafft --auto --thread 4 ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.fasta > ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.aln
mafft --auto --thread 4 ../trees/nifH/nifH_vnfH_anfH_clustered.fasta > ../trees/nifH/nifH_vnfH_anfH_clustered.aln

# remove gappy alignments
#trimal -in ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.aln -out ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -sgc -gappyout -keepheader
#trimal -in ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.aln -out ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.trim -sgc -gappyout -keepheader
trimal -in ../trees/nifH/nifH_vnfH_anfH_clustered.aln -out ../trees/nifH/nifH_vnfH_anfH_clustered.trim -sgc -gappyout -keepheader

echo "tree building"
# build maximum likelihood tree
#iqtree -s ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000 #use this to find best model and t>
#iqtree -s ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -safe -m LG+R10 -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000 #use this to find best model and threads
#iqtree -s ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000
iqtree -s ../trees/nifH/nifH_vnfH_anfH_clustered.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000

# Replace tree tip IDs with metadata-matched organism/cluster/genome/contig/operon IDs.
TREE_FILE="../trees/nifH/nifH_vnfH_anfH_clustered.trim.treefile"
python helper.py tree_node_match_metadata "$TREE_FILE"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
