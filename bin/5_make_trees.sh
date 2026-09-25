#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=18:15:00   # walltime # about 3hrs for ~300 seqs, 20+ hrs for 7000+
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 150GB   # memory per node
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
#cat ../results/final/fastas/final_nifH.fasta ../results/final/fastas/final_vnfH.fasta ../results/final/fastas/final_anfH.fasta > ../trees/nifH/nifH_anfH_vnfH.fasta
#cat ../results/final/fastas/final_nifD*.fasta ../results/final/fastas/final_anfD*.fasta ../results/final/fastas/final_vnfD*.fasta > ../diazoDB-comparison/tree-comparison/nifD_anfD_vnfD.fasta

#DIR="../diazoDB-comparison/tree-comparison"
GENE="H"
DIR="../trees/nif${GENE}"
TREE_FILE="${DIR}/nif${GENE}_anf${GENE}_vnf${GENE}.fasta"
CLUSTER="${DIR}/nif${GENE}_anf${GENE}_vnf${GENE}_clustered.fasta"

#mkdir -p "${DIR}/tmp"
#find "${DIR}/tmp"/ -type f -delete

#mmseqs createdb "$TREE_FILE" "${DIR}/tmp/seqDB"
#mmseqs cluster "${DIR}/tmp/seqDB" "${DIR}/tmp/clustered" "${DIR}/tmp" --min-seq-id 0.9 -c 0.8 --cov-mode 0
#mmseqs createtsv "${DIR}/tmp/seqDB" "${DIR}/tmp/seqDB" "${DIR}/tmp/clustered" "${CLUSTER}.tsv"
#mmseqs result2repseq "${DIR}/tmp/seqDB" "${DIR}/tmp/clustered" "${DIR}/tmp/clustered_reps"
#mmseqs result2flat "${DIR}/tmp/seqDB" "${DIR}/tmp/seqDB" "${DIR}/tmp/clustered_reps" "$CLUSTER" --use-fasta-header

# count clusters
#num=$(grep ">" "$CLUSTER" | wc -l)
#echo "$num clusters for 0.9"

# add outgroup
#cat ../trees/BchL.fasta >> "$CLUSTER"
#cat ../trees/CfbD.fasta ../trees/BchN.fasta ../trees/BchB.fasta >> "$CLUSTER"

# find comparison database closest match to clusters

#module load blast/2.15.0
#makeblastdb -in "$CLUSTER" -dbtype prot
#blastp -query ../diazoDB-comparison/Kacar-Results/swh:1:dir:7c3fb980f24a20df3144ebbf1f4be9feb80c634f/converted-nifD-extant.fasta -db "$CLUSTER" -out ../diazoDB-comparison/tree-comparison/NSDB.blast -outfmt 6 -max_target_seqs 5
#blastp -query ../diazoDB-comparison/NFixDB-Results/nifD_anfD_vnfD_12192023.faa -db "$CLUSTER" -out ../diazoDB-comparison/tree-comparison/NFixDB.blast -outfmt 6 -max_target_seqs 5
#blastp -query ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -db "$CLUSTER" -out ../diazoDB-comparison/tree-comparison/Nif-Finder.blast -outfmt 6 -max_target_seqs 5

# return sequences without a match in DiazoDB (pident > 90%)
#awk '$2 > 99 {print $1}' ../diazoDB-comparison/tree-comparison/NFixDB-full.blast | sort -u > ../diazoDB-comparison/tree-comparison/NFixDB-matches.txt
#seqkit grep -v -f ../diazoDB-comparison/tree-comparison/NFixDB-matches.txt ../diazoDB-comparison/NFixDB-Results/nifD_anfD_vnfD_12192023.faa -o ../diazoDB-comparison/tree-comparison/NFixDB-no_hits.fasta
#awk '$2 > 99 {print $1}' ../diazoDB-comparison/tree-comparison/NSDB-full.blast | sort -u > ../diazoDB-comparison/tree-comparison/NSDB-matches.txt
#seqkit grep -v -f ../diazoDB-comparison/tree-comparison/NSDB-matches.txt ../diazoDB-comparison/Kacar-Results/swh:1:dir:7c3fb980f24a20df3144ebbf1f4be9feb80c634f/converted-nifD-extant.fasta -o ../diazoDB-comparison/tree-comparison/NSDB-no_hits.fasta
#awk '$2 > 99 {print $1}' ../diazoDB-comparison/tree-comparison/Nif-Finder-full.blast | sort -u > ../diazoDB-comparison/tree-comparison/Nif-Finder-matches.txt
#seqkit grep -v -f ../diazoDB-comparison/tree-comparison/Nif-Finder-matches.txt ../diazoDB-comparison/Nif-finder-Results/true-nifD.faa -o ../diazoDB-comparison/tree-comparison/Nif-Finder-no_hits.fasta

#{
#awk '/^>/{print ">nifD|" substr($0,2); next} {print}' ../HMM_seeds/clustered_nifD_rep_seq.fasta
#awk '/^>/{print ">nifK|" substr($0,2); next} {print}' ../HMM_seeds/clustered_nifK_rep_seq.fasta
#awk '/^>/{print ">nifE|" substr($0,2); next} {print}' ../HMM_seeds/clustered_nifE_rep_seq.fasta
#awk '/^>/{print ">nifN|" substr($0,2); next} {print}' ../HMM_seeds/nifN_merged_len300.fasta
#awk '/^>/{print ">nifB|" substr($0,2); next} {print}' ../HMM_seeds/nifB.fasta
#} > ../HMM_seeds/nifDKENB.fasta

# align nif sequences
CLUSTER="${CLUSTER%.*}"
#mafft --auto --thread 4 ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.fasta > ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.aln
#mafft --auto --thread 4 ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.fasta > ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.aln
#mafft --auto --thread 4 ../trees/nifH/nifH_vnfH_anfH_clustered.fasta > ../trees/nifH/nifH_vnfH_anfH_clustered.aln
#mafft --auto --thread 4 "${CLUSTER}.fasta" > "${CLUSTER}.aln"

# remove gappy alignments
#trimal -in ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.aln -out ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -sgc -gappyout -keepheader
#trimal -in ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.aln -out ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.trim -sgc -gappyout -keepheader
#trimal -in ../trees/nifH/nifH_vnfH_anfH_clustered.aln -out ../trees/nifH/nifH_vnfH_anfH_clustered.trim -sgc -gappyout -keepheader
#trimal -in "${CLUSTER}.aln" -out "${CLUSTER}.trim" -sgc -gappyout -keepheader

echo "tree building"
# build maximum likelihood tree
#iqtree -s ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000 #use this to find best model and t>
#iqtree -s ../trees/nifK_noOut_04292025/clustered_nifK_noOut_rep_seq.trim -safe -m LG+R10 -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000 #use this to find best model and threads
#iqtree -s ../trees/nifH_500nodes/nifH_500nodes_clustered_rep_seq.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000
#iqtree -s ../trees/nifH/nifH_vnfH_anfH_clustered.trim -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000
iqtree -s "${CLUSTER}.trim" -pre "${CLUSTER}" -safe -m MFP -msub nuclear -T AUTO -ntmax 8 -B 1000 -alrt 1000

# Replace tree tip IDs with metadata-matched organism/cluster/genome/contig/operon IDs.
python helper.py tree_node_match_metadata "${CLUSTER}.treefile"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
