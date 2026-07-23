#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=05:12:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 200GB   # memory per CPU core
#SBATCH -J download   # job name

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# download GTDB representative protein sequences
curl -o protein_faa_reps_latest.tar.gz https://data.gtdb.aau.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.t
echo "rep genomes download done"

# download GTDB metadata files
curl https://data.gtdb.aau.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz
curl https://data.gtdb.aau.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz
echo "metadata download done"

# decompress the downloaded files
tar -xzf /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_latest.tar.gz
echo "tar done"

gunzip -r /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_latest/
gunzip -r /resnick/groups/enviromics/zahra/diazoDB-HPC/ar53_metadata.tsv.gz
gunzip -r /resnick/groups/enviromics/zahra/diazoDB-HPC/bac120_metadata.tsv.gz
echo "gunzip done"

# concatenate metadata files
cat /resnick/groups/enviromics/zahra/diazoDB-HPC/ar53_metadata.tsv /resnick/groups/enviromics/zahra/diazoDB-HPC/bac120_metadata.tsv > /resnick/groups/enviromics/zahra/diazoDB-HPC/bin/GTDB_metadata.tsv
gzip /resnick/groups/enviromics/zahra/diazoDB-HPC/bin/GTDB_metadata.tsv

# count clusters
num=$(ls -lh  /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_latest/archaea/*.faa | wc -l)
echo "$num files processed for archaea"

# count clusters
num=$(ls -lh  /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_latest/bacteria/ | wc -l)
echo "$num files processed for bacteria"

# add to calkit for remote storage + version control
calkit add /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_latest/* --to dvc-zip
echo "calkit done"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
