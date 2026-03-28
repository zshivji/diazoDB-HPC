#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=05:12:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 40GB   # memory per CPU core
#SBATCH -J download   # job name

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

#curl -O https://data.gtdb.aau.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz

#echo "download done"

tar -xzf /resnick/groups/enviromics/zahra/diazoDB-HPC/gtdb_proteins_aa_reps.tar.gz

echo "tar done"

gunzip -r /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps/

echo "gunzip done"

# count clusters
num=$(ls -lh  /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps/archaea/*.faa | wc -l)
echo "$num files processed for archaea"

# count clusters
num=$(ls -lh  /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps/bacteria/ | wc -l)
echo "$num files processed for bacteria"

calkit add /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps/* -t dvc

echo "calkit done"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
