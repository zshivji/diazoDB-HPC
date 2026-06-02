import os, base64, subprocess, requests, logging
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

API       = os.environ["API_URL"].rstrip("/")
HEADERS   = {"x-runner-token": os.environ["RUNNER_SECRET"]}
JOB_BASE  = Path(os.environ["HPC_JOB_BASE"])
DUMMY_RUNNER = os.environ.get("DUMMY_RUNNER", "").lower() in {"1", "true", "yes"}
SLURM_SCRIPT = os.environ.get("SLURM_SCRIPT")
SLURM_LOG_DIR = os.environ.get("SLURM_LOG_DIR")

# runner polls API (HTTP request)
def poll() -> list[dict]:
    r = requests.get(f"{API}/api/v1/runner/jobs", headers=HEADERS, timeout=15)
    print(r.status_code)
    print(r.text)
    r.raise_for_status()
    return r.json()

# marks job as processing/failed/complete
def mark(job_id: str, status: str, error: str = "") -> None:
    payload = {"status": status}
    if error:
        payload["error_message"] = error
    requests.patch(
        f"{API}/api/v1/runner/jobs/{job_id}",
        json=payload,
        headers=HEADERS,
        timeout=10,
    )

# uploads analysis output back to server
def push_result(job_id: str, result_path: Path) -> None:
    content_types = {
        ".csv": "text/csv",
        ".html": "text/html",
        ".htm": "text/html",
        ".pdf": "application/pdf",
    }
    content_type = content_types.get(result_path.suffix.lower(), "application/octet-stream")
    data = base64.b64encode(result_path.read_bytes()).decode()
    r = requests.post(
        f"{API}/api/v1/runner/jobs/{job_id}/result",
        json={"filename": result_path.name, "content_type": content_type, "data_base64": data},
        headers=HEADERS,
        timeout=60,
    )
    r.raise_for_status()


def download_input(job: dict, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(
        f"{API}/api/v1/runner/jobs/{job['id']}/input",
        headers=HEADERS,
        timeout=120,
    )
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


def workspace_for(job_id: str) -> Path:
    return JOB_BASE / job_id


# def resolve_input_path(job: dict) -> Path:
#     job_id = job["id"]
#     filename = job["filename"]
#     preferred = workspace_for(job_id) / filename
#     if preferred.exists():
#         return preferred

#     hpc_path = Path(job["hpc_path"])
#     if hpc_path.exists():
#         return hpc_path

#     if hpc_path.is_absolute() and JOB_BASE.name in hpc_path.parts:
#         base_index = hpc_path.parts.index(JOB_BASE.name)
#         relative_parts = hpc_path.parts[base_index + 1:]
#         return JOB_BASE.joinpath(*relative_parts)

#     return preferred


def run_slurm(job_id: str, input_path: Path, out_file: Path) -> None:
    if not SLURM_SCRIPT:
        raise RuntimeError("SLURM_SCRIPT is not set")

    out_file.parent.mkdir(parents=True, exist_ok=True)
    if SLURM_LOG_DIR:
        Path(SLURM_LOG_DIR).mkdir(parents=True, exist_ok=True)
    export_vars = ",".join([
        "ALL",
        f"DIAZODB_JOB_ID={job_id}",
        f"DIAZODB_INPUT={input_path}",
        f"DIAZODB_OUTDIR={out_file.parent}",
        f"DIAZODB_OUTPUT={out_file}",
        f"DIAZODB_USE_PRODIGAL={os.environ.get('DIAZODB_USE_PRODIGAL', 'false')}",
    ])
    subprocess.run(
        [
            "sbatch",
            "--wait",
            "--parsable",
            "--job-name",
            f"diazodb_{job_id[:8]}",
            "--export",
            export_vars,
            SLURM_SCRIPT,
        ],
        check=True,
        timeout=int(os.environ.get("SLURM_WAIT_TIMEOUT_SECONDS", "7200")), # currently set at 2hrs
    )


def run_analysis(job_id: str, input_path: Path, out_file: Path) -> None:
    if DUMMY_RUNNER: # for testing
        out_file.write_text(
            "sequence_id,prediction,confidence\n"
            f"{input_path.name},dummy_result,1.0\n",
            encoding="utf-8",
        )
        return

    if SLURM_SCRIPT:
        run_slurm(job_id, input_path, out_file)
        return

    # # ── Replace this command with your actual non-Slurm analysis tool ──
    # subprocess.run(
    #     ["your_tool", str(input_path), "--output", str(out_file)],
    #     check=True,
    #     timeout=7200,
    # )


def process(job: dict) -> None:
    job_id = job["id"]
    os.environ["DIAZODB_USE_PRODIGAL"] = str(job.get("use_prodigal", False)).lower()
    # input_path = resolve_input_path(job)
    # out_file = input_path.parent / "output.csv"

    input_path = download_input(job, workspace_for(job_id) / job["filename"])
    out_file = input_path.parent / "output.csv"

    if not input_path.exists():
        hpc_path = Path(job["hpc_path"])
        log.error(f"[{job_id}] input not found: {hpc_path} (resolved to {input_path})")
        mark(job_id, "failed", "Input file not found on HPC")
        return

    mark(job_id, "processing") # update status
    log.info(f"[{job_id}] processing {input_path}")

    try:
        run_analysis(job_id, input_path, out_file)
        if not out_file.exists():
            raise FileNotFoundError(f"Expected result file was not created: {out_file}")
        push_result(job_id, out_file)
        log.info(f"[{job_id}] complete")
    except subprocess.CalledProcessError as e:
        mark(job_id, "failed", str(e))
        log.error(f"[{job_id}] tool failed: {e}")
    except Exception as e:
        mark(job_id, "failed", str(e))
        log.error(f"[{job_id}] unexpected error: {e}")


if __name__ == '__main__':
    log.info(f"Runner started")
    try:
        jobs = poll()
        log.info(f"Polled — {len(jobs)} job(s) ready")
        for job in jobs:
            process(job)
    except Exception as e:
        log.error(f"Poll error: {e}")
