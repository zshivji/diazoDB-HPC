DATASETS

NFixDB (updated April 9)
- published: June 6, 2024
- paper: https://doi.org/10.1093/nargab/lqae063
- data: https://zenodo.org/records/15177936 (fasta files dated nifX_12192023.faa)
- License: 
- notes: 
  - do not group, but differentiate between nif/anf/vnf
  - no nifE, nifN, nifB
  - also searched in GTDB
  - filtering criteria
    - aligment length > 150 a.a.
    - e-value < 9.9 × 10−15

- number of sequences
    nifH: 4085
    nifD: 4162
    nifK: 4158
    anfH: 261
    anfD: 263
    anfK: 262
    vnfH: 56
    vnfD: 59
    vnfK: 59
- groups
    nif: 4085
    anf: 261
    vnf: 56

NSDB (updated April 9)
- published: Sept 11, 2025
- paper: https://elifesciences.org/articles/105613#content
- website: https://nsdb.bact.wisc.edu
- data: https://archive.softwareheritage.org/browse/origin/directory/?origin_url=https://github.com/kacarlab/nif-structures
- License: MIT
- notes:
    - only nifH , nifD , nifK alignment files available
    - no nifE, nifN, nifB (?)
    - full results not avail on website unless search (search "_" to view Anc and extant)
    - total DB entries is 5,378
        - but this includes separate entries for DDKK and HH
        - and 6x alternates for each ancestral sequence
        - (385 unique extant + 6*384 unique ancestral ) * 2 = 5378
    - can you download tabular metadata?
        - on data archive in data/tree/ have nif groups in clade-asignments.csv for 764 organisms (5 entries missing from ancestral)
- number of sequences
    extant nif: 385
- groups
    nif: 351
    (nif-i: 139)
    (nif-ii: 181)
    (nif-iii: 31)
    anf: 17
    vnf: 17

Cyano/Nif-Finder (updated April 9)
- published: BioRxiv, Jan 15 2026
- paper: https://doi.org/10.64898/2026.01.15.699626
- data:
- License: MIT
- notes:
    - nifHDKENB
- number of sequences
    - 285 nifHDKENB
- groups
    nif: 351
    anf: 
    vnf: 

TOOLS

NifFinder
- published: Oct 2025
- paper: https://doi.org/10.1093/bioadv/vbaf260
- data: N/A
_ license: 
- notes:
    - software tool not dataset
    - how do they define accuracy & sensitivity? what is the test set
- number of sequences
    nifH:
    nifD:
    nifK:
- groups: N/A


Carmna:
- published: Mar 2025
- paper: https://doi.org/10.1093/bib/bbaf197
- data: N/A
_ license: N/A
- notes:
    - software tool not dataset
    - how do they define accuracy & sensitivity? what is the test set
- number of sequences
    nifH:
    nifD:
    nifK:
- groups: N/A

NifPred:
- published: May 2018
- paper: https://doi.org/10.1093/bib/bbaf197
- data: 
_ license: 
- notes:
    - reviewed by former orphan lab member
- number of sequences
    nifH:
    nifD:
    nifK:
    nifE:
    nifN:
    nifB:
- groups: N/A