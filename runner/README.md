# DiazoDB HPC Runner

The runner is intended to run on the HPC login/service node. It polls the web
API for ready jobs, downloads each submitted input file over HTTPS, submits the
analysis to Slurm, and posts the final result back to the API.

This avoids SSH/SCP from the web server into the cluster. The only secret shared
between systems is `RUNNER_SECRET`. Set `RUNNER_PULL_ONLY=true` in the web
API .env to skip Globus transfers and let the runner download inputs.

## Minimal Environment

```bash
API_URL=https://api.example.edu
RUNNER_SECRET=<same value as website RUNNER_SECRET>
HPC_JOB_BASE=/scratch/zshivji/diazodb/jobs
SLURM_SCRIPT=/path/to/diazoDB-HPC/bin/diazodb_classify.sh
SLURM_WAIT_TIMEOUT_SECONDS=32400
```

Optional pipeline settings consumed by `bin/diazodb_classify.sh`:

```bash
DIAZODB_CONDA_ENV=/path/to/conda/env
DIAZODB_HMM_PROFILE=/path/to/combined_nif.hmm
DIAZODB_USE_PRODIGAL=false
```

For smoke testing without Slurm:

```bash
DUMMY_RUNNER=true python runner.py
```

The runner performs one poll per invocation and is intended to be invoked
periodically by cron.

## Per-job files

Each API job is isolated under its UUID:

```text
HPC_JOB_BASE/<job-id>/
├── input/                         # downloaded user input
├── intermediate/
│   ├── predicted_proteins.faa    # present when Prodigal is requested
│   ├── hmmsearch_results/
│   │   ├── hmm_out/               # raw .out and .domtblout files
│   │   └── hits.csv               # this job's parsed HMM hits
│   ├── conserved_res/             # FASTAs, alignments, residue checks
│   └── proteins/
│       ├── user/                   # uploaded or Prodigal-predicted proteins
│       └── reference/              # packaged conserved-residue references
├── results/
│   ├── nif_clusters.csv           # file uploaded to the API
│   ├── nif_final.csv              # per-gene conserved-residue calls
│   └── fastas/                    # final FASTAs grouped by gene call
└── logs/
    ├── slurm-<name>-<id>.out
    └── slurm-<name>-<id>.err
```

The runner passes the job-specific parsed hits, temporary-results directory,
final-results directory, protein lookup directory, and reference IDs directly
to `conserved-res.py`. Separate cron invocations therefore do not share mutable
analysis files.

Runner jobs do not read the GTDB representative-genome corpus or GTDB metadata.
They use the uploaded genome/proteome plus the small fixed reference FASTA in
`runner/references/`. Bacterial and archaeal uploads follow the same execution
path.

Input downloads are streamed to disk in 8 MiB chunks, so the runner does not
hold a multi-gigabyte upload in memory.
