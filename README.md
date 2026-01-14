# DiazoDB

## Nitrogenase Database
Automated annotation and curation of nitrogenase genes using profile hidden Markov models and conserved residue matching.

## Scripts

+ **1_align_mafft.sh** - Clusters (mmseqs easy-cluster) HMM seed sequences at 90% AAI. Aligns (MAFFT) clustered seeds sequences. Builds HMM profiles and combines into single file.

+ **2_hmmsearch.sh** - Runs pHMM search against GTDB R220 all_rep_proteins_aa database (bacteria and archaea).

+ **3_parse_hmm.sh**

+ **4_conserved_res.sh**

+ **5_make_trees.sh**

+ **6_operon_org.sh**

+ **7_SSN.sh**

+ **Parse_hmm_results.py**

+ **Parse_tophits.py**

+ **aln_nif_hits.py**

+ **cluster_pos.py**

+ **conserved-res.py**

+ **final-fasta-export.py**

+ **get-operon.py**

+ **operon-org-plot.py**

+ **run_check_seq.sh**



## examples

+ **final.py** - makes/formats final TSV, outputs ***NFixDB.tsv***

+ **nitrogenase_fastas.py** - creates new fasta files from ***filteredfasta.tsv***

