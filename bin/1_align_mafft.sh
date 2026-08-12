#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=00:12:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH --job-name=hmm_seed   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

module load mafft/7.505-gcc-13.2.0-nklkvtc

# cluster
mmseqs easy-cluster ../HMM_seeds/nifH_merged_len200.fasta ../HMM_seeds/clustered_nifH tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/nifD_merged_len300.fasta ../HMM_seeds/clustered_nifD tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/nifK_merged_len300.fasta ../HMM_seeds/clustered_nifK tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/nifE_merged_len300.fasta ../HMM_seeds/clustered_nifE tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/nifN_merged_len300.fasta ../HMM_seeds/clustered_nifN tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/nifB.fasta ../HMM_seeds/clustered_nifB tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/anfO.fasta ../HMM_seeds/clustered_anfO tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/anfG.fasta ../HMM_seeds/clustered_anfG tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
mmseqs easy-cluster ../HMM_seeds/vnfG.fasta ../HMM_seeds/clustered_vnfG tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0

# align nif sequences
mafft --auto --thread 4 ../HMM_seeds/clustered_nifH_rep_seq.fasta > ../alignments/nifH-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_nifD_rep_seq.fasta > ../alignments/nifD-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_nifK_rep_seq.fasta > ../alignments/nifK-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_nifE_rep_seq.fasta > ../alignments/nifE-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_nifN_rep_seq.fasta > ../alignments/nifN-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_nifB_rep_seq.fasta > ../alignments/nifB-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_anfO_rep_seq.fasta > ../alignments/anfO-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_anfG_rep_seq.fasta > ../alignments/anfG-07292026.faa
mafft --auto --thread 4 ../HMM_seeds/clustered_vnfG_rep_seq.fasta > ../alignments/vnfG-07292026.faa

# rm exisiting profiles
rm ./HMMs/*

# build hmm profiles
hmmbuild ../HMMs/nifH_07292026.hmm ../alignments/nifH-07292026.faa
hmmbuild ../HMMs/nifD_07292026.hmm ../alignments/nifD-07292026.faa
hmmbuild ../HMMs/nifK_07292026.hmm ../alignments/nifK-07292026.faa
hmmbuild ../HMMs/nifE_07292026.hmm ../alignments/nifE-07292026.faa
hmmbuild ../HMMs/nifN_07292026.hmm ../alignments/nifN-07292026.faa
hmmbuild ../HMMs/nifB_07292026.hmm ../alignments/nifB-07292026.faa
hmmbuild ../HMMs/anfO_07292026.hmm ../alignments/anfO-07292026.faa
hmmbuild ../HMMs/anfG_07292026.hmm ../alignments/anfG-07292026.faa
hmmbuild ../HMMs/vnfG_07292026.hmm ../alignments/vnfG-07292026.faa

# combine hmm profiles
cat ../HMMs/*.hmm > ../HMMs/combined_nif_07292026.hmm

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
