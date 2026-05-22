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
SLURM_SCRIPT=/path/to/diazoDB-HPC/runner/slurm/diazodb_classify.sbatch
SLURM_WAIT_TIMEOUT_SECONDS=7200
```

Optional pipeline settings consumed by `runner/slurm/diazodb_classify_pipeline.sh`:

```bash
DIAZODB_CONDA_ENV=/path/to/conda/env
DIAZODB_HMM_PROFILE=/path/to/combined_nif.hmm
DIAZODB_USE_PRODIGAL=false
DIAZODB_PARSE_COMMAND='python /path/to/parse_uploaded_job.py --input "$QUERY_FASTA" --tblout "$OUTDIR/hmmsearch.tbl" --output "$FINAL_OUTPUT"'
```

For smoke testing without Slurm:

```bash
DUMMY_RUNNER=true python runner.py
```

For smoke testing with Slurm but placeholder parser output, set `SLURM_SCRIPT`
and leave `DIAZODB_PARSE_COMMAND` unset.
