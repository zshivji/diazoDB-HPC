#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=00:12:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH -J mafft.e%j   # job name

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

module load mafft/7.505-gcc-13.2.0-nklkvtc

# cluster
#mmseqs easy-cluster ../input-files/nifH_merged_len200.fasta ../input-files/clustered_nifH tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
#mmseqs easy-cluster ../input-files/nifD_merged_len300.fasta ../input-files/clustered_nifD tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
#mmseqs easy-cluster ../input-files/nifK_merged_len300.fasta ../input-files/clustered_nifK tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
#mmseqs easy-cluster ../input-files/nifE_merged_len300.fasta ../input-files/clustered_nifE tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
#mmseqs easy-cluster ../input-files/nifN_merged_len300.fasta ../input-files/clustered_nifN tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
#mmseqs easy-cluster ../input-files/nifB_merged_len100.fasta ../input-files/clustered_nifB tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
#mmseqs easy-cluster ../input-files/nifV-pfam-blast.fasta ../input-files/clustered_nifV-pfam-blast tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0

# align nif sequences
#mafft --auto --thread 4 ../input-files/clustered_nifH_rep_seq.fasta > ../alignments/nifH-03292025.faa
#mafft --auto --thread 4 ../input-files/clustered_nifD_rep_seq.fasta > ../alignments/nifD-03292025.faa
#mafft --auto --thread 4 ../input-files/clustered_nifK_rep_seq.fasta > ../alignments/nifK-03292025.faa
#mafft --auto --thread 4 ../input-files/clustered_nifE_rep_seq.fasta > ../alignments/nifE-03292025.faa
#mafft --auto --thread 4 ../input-files/clustered_nifN_rep_seq.fasta > ../alignments/nifN-03292025.faa
#mafft --auto --thread 4 ../input-files/clustered_nifB_rep_seq.fasta > ../alignments/nifB-03292025.faa
#mafft --auto --thread 4 ../input-files/clustered_nifV-pfam-blast_rep_seq.fasta > ../alignments/nifV-07242025.faa

eval "$(conda shell.bash hook)"
conda activate /central/groups/enviromics/miniconda3/envs/hmmer

# build hmm profiles
#hmmbuild ../HMMs/nifH_03292025.hmm ../alignments/nifH-03292025.faa
#hmmbuild ../HMMs/nifD_03292025.hmm ../alignments/nifD-03292025.faa
#hmmbuild ../HMMs/nifK_03292025.hmm ../alignments/nifK-03292025.faa
#hmmbuild ../HMMs/nifE_03292025.hmm ../alignments/nifE-03292025.faa
#hmmbuild ../HMMs/nifN_03292025.hmm ../alignments/nifN-03292025.faa
#hmmbuild ../HMMs/nifB_03292025.hmm ../alignments/nifB-03292025.faa
#hmmbuild ../HMMs/nifV_07242025.hmm ../alignments/nifV-07242025.faa
#hmmbuild ../HMMs/pchlide_04122025.hmm ../alignments/pchlide-04012025.faa

# combine hmm profiles
#cat ../HMMs/*.hmm > ../HMMs/combined_nif_04012025.hmm

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
