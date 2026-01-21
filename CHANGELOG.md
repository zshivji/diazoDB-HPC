# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### To-do
- update nitrogenase filtering criteria
- differentiate anf, vnf from nif


### Added

- Description of scripts in repo to README.md
- Sandbox diazoDB website files under /website
- variable **date** so all files in parse_hmmer pipeline automatically named by date created
- need test case for keeping/dropping best hit in each gene cluster for conserved residue matching
 
### Fixed 
- compress GTDB-taxonomy.tsv and update code to properly read in
- Major upgrades to parse_hmmer.ipynb
    - changed # simultaneously aligned seqs from 4000 to 200

    - modified nifD, nifE conserved residue matching
    - removed drop duplicates (without specifying same hit) from conserved residue matching
- Updated aln_nif_hits.py to reflect parse_hmmer.ipynb changes

### Changed

### Removed


## [0.0.1] - 2026-01-12

### Inital Commit

- Uploaded code, README.md, and CHANGELOG.md to github repo https://github.com/zshivji/diazoDB-HPC
