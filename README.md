# DiazoDB

## Nitrogenase Database
Automated annotation and curation of nitrogenase genes using profile hidden Markov models and conserved residue matching.

---

## Repo structure

```
diazo-db/
│
├── index.html                  ← Homepage (article listing)
├── about.html                  ← About page (create when ready)
├── residues.html               ← Residue index page (create when ready)
│
├── css/
│   ├── main.css                ← Shared styles: header, footer, homepage layout
│   └── article.css             ← Article-page styles: TOC, body typography, tables
│
├── js/
│   ├── articles.js             ← Article registry + homepage renderer
│   ├── article.js              ← TOC scroll-spy (loaded only on article pages)
│   └── main.js                 ← Search and tab filtering (loaded on homepage)
│
├── images/
│   ├── thumbs/                 ← Homepage card thumbnails (88×66 px recommended)
│   │   └── example.jpg
│   └── figures/                ← Full-size figures used inside articles
│       └── example.png
│
└── articles/
    ├── article-template/       ← COPY THIS to start a new article
    │   └── index.html
    │
    ├── femo-cofactor/
    │   └── index.html
    ├── mofe-protein/
    │   └── index.html
    ├── nifa-sigma54/
    │   └── index.html
    └── ...
```

---

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

---

## Typography

| Use | Font |
|-----|------|
| Body text, UI, nav | DM Sans |
| Article titles, lead paragraphs | Lora (serif) |
| Residue names, code | DM Mono |

Both fonts load from Google Fonts. If you want self-hosted fonts for offline use,
download from [fonts.google.com](https://fonts.google.com) and update the `@font-face`
rules in `css/main.css`.