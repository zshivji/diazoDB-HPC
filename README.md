# DiazoDB

## Nitrogenase Database
Automated annotation and curation of nitrogenase genes using profile hidden Markov models and conserved residue matching.

## Install

## Scripts

+ **1_align_mafft.sh** - Clusters (mmseqs easy-cluster) HMM seed sequences at 90% AAI. Aligns (MAFFT) clustered seeds sequences. Builds HMM profiles and combines into single file. Stores alignmnets in ***alignments/*** and profiles in ***HMMs/***.

+ **2_hmmsearch.sh** - Runs pHMM search against GTDB R220 all_rep_proteins_aa database (bacteria and archaea). Stores results as individual files in ***results/{archaea | bacteria}*** (per GTDB accession number) 

+ **3_parse_hmm.sh** - Calls **parse_hmm.py** to parse HMM hits and filter top hits.

+ **4_conserved_res.sh** - Calls (1) ***aln_nif_hits.py***, (2) ***conserved-res.py***, (3) ***final-fasta-export.py***, and (4) ***diazoDB-check.py***.

+ **5_make_trees.sh** - 

+ **6_operon_org.sh**

+ **7_SSN.sh**

+ **parse_hmm.py** - Parse pHMM search results (see below criteria), combines with GTDB taxonomy data, stores results in ***hits.feather***, and can filter top hits in one call.
    + positive bit score
    + full sequence evalue must be significant (<0.01)
    + best domain evalue should be significant (<0.01)
        + otherwise flagged for manual review to check if it is distant homolog or just short repeats

+ **aln_nif_hits.py**

+ **cluster_pos.py**

+ **conserved-res.py**

+ **final-fasta-export.py**

+ **diazoDB-check.sh**

## Data Provenance and Reproducibility

### 1. Source Data
**Source Data:** Genome Taxonomy Database (GTDB) R232  
**Source Files:**
+ https://data.gtdb.aau.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz  
+ https://data.gtdb.aau.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz
+ https://data.gtdb.aau.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz

### 2. Manual Curation 
The following files were manually curated:
+ HMM_seeds/nifH.fasta
+ HMM_seeds/nifD.fasta
+ HMM_seeds/nifK.fasta
+ HMM_seeds/nifE.fasta
+ HMM_seeds/nifN.fasta
+ HMM_seeds/nifB.fasta
+ HMM_seeds/anfG.fasta
+ HMM_seeds/anfO.fasta
+ HMM_seeds/vnfG.fasta

+ **get-operon.py**

+ **operon-org-plot.py**


