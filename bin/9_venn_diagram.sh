#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=2:15:00   # walltime # about 3hrs for ~300 seqs, 20+ hrs for 7000+
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 10GB   # memory per node
#SBATCH --job-name=venn   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# load env, software
eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/seq-tools

# NIFD, to keep full fasta header, run easy-NIFD workflow separately
#cat ../results/final/fastas/final_nifD*.fasta ../results/final/fastas/final_anfD*.fasta ../results/final/fastas/final_vnfD*.fasta > "${VENN_DIR}nifD_anfD_vnfD.fasta
NIFD="../diazoDB-comparison/tree-comparison/nifD_anfD_vnfD.fasta"
VENN_DIR="../diazoDB-comparison/venn"
DB="${VENN_DIR}/nifD_anfD_vnfD.fasta"

# find comparison database closest match to NIFDs
module load blast/2.15.0
#makeblastdb -in "$NIFD" -dbtype prot -out "$DB" 
#blastp -query ../diazoDB-comparison/Kacar-Results/swh:1:dir:7c3fb980f24a20df3144ebbf1f4be9feb80c634f/converted-nifD-extant.fasta -db "$DB" -out "${VENN_DIR}/NSDB.blast" -outfmt 6 -max_target_seqs 5
#blastp -query ../diazoDB-comparison/NFixDB-Results/nifD_anfD_vnfD_12192023.faa -db "$DB" -out "${VENN_DIR}/NFixDB.blast" -outfmt 6 -max_target_seqs 5
#blastp -query ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -db "$DB" -out "${VENN_DIR}/Nif-Finder.blast" -outfmt 6 -max_target_seqs 5

# return sequences without a match in DiazoDB (pident > 97.5%)
awk '$3 > 97.5 {print $1}' "${VENN_DIR}/NFixDB.blast" | sort -u > "${VENN_DIR}/NFixDB-matches.txt"
seqkit grep -v -f "${VENN_DIR}/NFixDB-matches.txt" ../diazoDB-comparison/NFixDB-Results/nifD_anfD_vnfD_12192023.faa -o "${VENN_DIR}/NFixDB-no_hits.fasta"
awk '$3 > 97.5 {print $1}' "${VENN_DIR}/NSDB.blast" | sort -u > "${VENN_DIR}/NSDB-matches.txt"
seqkit grep -v -f "${VENN_DIR}/NSDB-matches.txt" ../diazoDB-comparison/Kacar-Results/swh:1:dir:7c3fb980f24a20df3144ebbf1f4be9feb80c634f/converted-nifD-extant.fasta -o "${VENN_DIR}/NSDB-no_hits.fasta"
awk '$3 > 97.5 {print $1}' "${VENN_DIR}/Nif-Finder.blast" | sort -u > "${VENN_DIR}/Nif-Finder-matches.txt"
seqkit grep -v -f "${VENN_DIR}/Nif-Finder-matches.txt" ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -o "${VENN_DIR}/Nif-Finder-no_hits.fasta"

#inter-database comparison
DB="${VENN_DIR}/NSDB"
#makeblastdb -in ../diazoDB-comparison/Kacar-Results/swh:1:dir:7c3fb980f24a20df3144ebbf1f4be9feb80c634f/converted-nifD-extant.fasta -dbtype prot -out "$DB"
#blastp -query ../diazoDB-comparison/NFixDB-Results/nifD_anfD_vnfD_12192023.faa -db "$DB" -out "${VENN_DIR}/NSDB-v-NFixDB.blast" -outfmt 6 -max_target_seqs 5
#blastp -query ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -db "$DB" -out "${VENN_DIR}/NSDB-v-Nif-Finder.blast" -outfmt 6 -max_target_seqs 5
awk '$3 > 97.5 {print $1}' "${VENN_DIR}/NSDB-v-NFixDB.blast" | sort -u > "${VENN_DIR}/NSDB-v-NFixDB-matches.txt"
seqkit grep -v -f "${VENN_DIR}/NSDB-v-NFixDB-matches.txt" ../diazoDB-comparison/NFixDB-Results/nifD_anfD_vnfD_12192023.faa -o "${VENN_DIR}/NSDB-v-NFixDB-no_hits.fasta"
awk '$3 > 97.5 {print $1}' "${VENN_DIR}/NSDB-v-Nif-Finder.blast" | sort -u > "${VENN_DIR}/NSDB-v-Nif-Finder-matches.txt"
seqkit grep -v -f "${VENN_DIR}/NSDB-v-Nif-Finder-matches.txt" ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -o "${VENN_DIR}/NSDB-v-Nif-Finder-no_hits.fasta"

DB="${VENN_DIR}/NFixDB"
#makeblastdb -in ../diazoDB-comparison/NFixDB-Results/nifD_anfD_vnfD_12192023.faa -dbtype prot -out "$DB"
#blastp -query ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -db "$DB" -out "${VENN_DIR}/NFixDB-v-Nif-Finder.blast" -outfmt 6 -max_target_seqs 5
awk '$3 > 97.5 {print $1}' "${VENN_DIR}/NFixDB-v-Nif-Finder.blast" | sort -u > "${VENN_DIR}/NFixDB-v-Nif-Finder-matches.txt"
seqkit grep -v -f "${VENN_DIR}/NFixDB-v-Nif-Finder-matches.txt" ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -o "${VENN_DIR}/NFixDB-v-Nif-Finder-no_hits.fasta"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
