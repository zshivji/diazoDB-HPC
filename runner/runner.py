import os, base64, json, shutil, subprocess, requests, logging
from datetime import datetime, timezone
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
HPC_DATABASE_BASE = os.environ.get("HPC_DATABASE_BASE")

# runner polls API (HTTP request)
def poll() -> list[dict]:
    r = requests.get(f"{API}/api/v1/runner/jobs", headers=HEADERS, timeout=15)
    log.info("API response status: %s", r.status_code)
    log.info("API response body: %s", r.text)
    r.raise_for_status()
    return r.json()

# marks job as processing/failed/complete
def mark(job_id: str, status: str, error: str = "") -> dict:
    payload = {"status": status}
    if error:
        payload["error_message"] = error
    r = requests.patch(
        f"{API}/api/v1/runner/jobs/{job_id}",
        json=payload,
        headers=HEADERS,
        timeout=10,
    )
    r.raise_for_status()
    response = r.json()
    if status == "failed" and not response.get("email_sent"):
        log.error(
            "[%s] job marked failed, but failure email was not sent: %s",
            job_id,
            response.get("email_status", "unknown"),
        )
    return response

# uploads analysis output back to server
def push_results(job_id: str, result_paths: list[Path]) -> dict:
    content_types = {
        ".csv": "text/csv",
        ".fa": "text/plain",
        ".faa": "text/plain",
        ".fasta": "text/plain",
        ".fna": "text/plain",
        ".html": "text/html",
        ".htm": "text/html",
        ".pdf": "application/pdf",
        ".json": "application/json",
        ".png": "image/png",
        ".svg": "image/svg+xml",
    }
    files = []
    for result_path in result_paths:
        content_type = content_types.get(
            result_path.suffix.lower(), "application/octet-stream"
        )
        files.append({
            "filename": result_path.name,
            "content_type": content_type,
            "data_base64": base64.b64encode(result_path.read_bytes()).decode(),
        })
    r = requests.post(
        f"{API}/api/v1/runner/jobs/{job_id}/result",
        json={"files": files},
        headers=HEADERS,
        timeout=60,
    )
    r.raise_for_status()
    response = r.json()
    if not response.get("email_sent"):
        log.error(
            "[%s] result stored by API, but notification email was not sent: %s",
            job_id,
            response.get("email_status", "unknown"),
        )
    return response


def push_result(job_id: str, result_path: Path) -> dict:
    """Backward-compatible helper for callers that have one result file."""
    return push_results(job_id, [result_path])


def archive_database_job(
    job: dict, input_path: Path, result_paths: list[Path], workspace: Path
) -> Path | None:
    """Persist opted-in inputs and results outside temporary HPC storage."""
    if not job.get("include_in_database"):
        return None
    if not HPC_DATABASE_BASE:
        raise RuntimeError("HPC_DATABASE_BASE is required for database contributions")

    archive_root = Path(HPC_DATABASE_BASE)
    archive_root.mkdir(parents=True, exist_ok=True)
    archive_path = archive_root / str(job["id"])
    temporary_path = archive_root / f".{job['id']}.tmp"
    if archive_path.exists():
        raise FileExistsError(f"Database archive already exists: {archive_path}")
    if temporary_path.exists():
        shutil.rmtree(temporary_path)

    try:
        (temporary_path / "input").mkdir(parents=True)
        (temporary_path / "results").mkdir(parents=True)
        shutil.copy2(input_path, temporary_path / "input" / input_path.name)

        archived_results = []
        for result_path in result_paths:
            relative_path = result_path.relative_to(workspace)
            destination = temporary_path / relative_path
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(result_path, destination)
            archived_results.append(str(relative_path))

        manifest = {
            "job_id": str(job["id"]),
            "filename": job["filename"],
            "orcid": job.get("orcid"),
            "archived_at": datetime.now(timezone.utc).isoformat(),
            "input": str(Path("input") / input_path.name),
            "results": archived_results,
        }
        (temporary_path / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        temporary_path.rename(archive_path)
    except Exception:
        shutil.rmtree(temporary_path, ignore_errors=True)
        raise

    return archive_path


def download_input(job: dict, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(
        f"{API}/api/v1/runner/jobs/{job['id']}/input",
        headers=HEADERS,
        stream=True,
        timeout=(15, 3600),
    ) as r:
        r.raise_for_status()
        with dest.open("wb") as output:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    output.write(chunk)
    return dest


def workspace_for(job_id: str) -> Path:
    return JOB_BASE / job_id


def run_slurm(
    job_id: str,
    input_path: Path,
    intermediate_dir: Path,
    result_path: Path,
    log_dir: Path,
) -> None:
    if not SLURM_SCRIPT:
        raise RuntimeError("SLURM_SCRIPT is not set")

    intermediate_dir.mkdir(parents=True, exist_ok=True)
    result_path.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    export_vars = ",".join([
        "ALL",
        f"DIAZODB_JOB_ID={job_id}",
        f"DIAZODB_INPUT={input_path}",
        f"DIAZODB_OUTDIR={intermediate_dir}",
        f"DIAZODB_OUTPUT={result_path}",
        f"DIAZODB_USE_PRODIGAL={os.environ.get('DIAZODB_USE_PRODIGAL', 'false')}",
    ])
    subprocess.run(
        [
            "sbatch",
            "--wait",
            "--parsable",
            "--job-name",
            f"diazodb_{job_id[:8]}",
            "--output",
            str(log_dir / "%x-%j.out"),
            "--export",
            export_vars,
            SLURM_SCRIPT,
        ],
        check=True,
        timeout=int(os.environ.get("SLURM_WAIT_TIMEOUT_SECONDS", "32400")),
    )


def run_analysis(
    job_id: str,
    input_path: Path,
    intermediate_dir: Path,
    result_path: Path,
    log_dir: Path,
) -> None:
    if DUMMY_RUNNER: # for testing
        result_path.parent.mkdir(parents=True, exist_ok=True)
        result_path.write_text(
            "sequence_id,prediction,confidence\n"
            f"{input_path.name},dummy_result,1.0\n",
            encoding="utf-8",
        )
        return

    if SLURM_SCRIPT:
        run_slurm(job_id, input_path, intermediate_dir, result_path, log_dir)
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

    workspace = workspace_for(job_id)
    safe_filename = Path(job["filename"].replace("\\", "/")).name
    if not safe_filename or safe_filename in {".", ".."}:
        raise ValueError("Job filename is invalid")
    input_path = download_input(job, workspace / "input" / safe_filename)
    intermediate_dir = workspace / "intermediate"
    result_path = workspace / "results" / "nif_clusters.csv"
    final_result_path = workspace / "results" / "nif_final.csv"
    operon_plot_path = workspace / "intermediate" / "operon-org.png"
    operon_metadata_path = workspace / "results" / "operon_metadata.json"
    log_dir = workspace / "logs"

    if not input_path.exists():
        hpc_path = Path(job["hpc_path"])
        log.error(f"[{job_id}] input not found: {hpc_path} (resolved to {input_path})")
        mark(job_id, "failed", "Input file not found on HPC")
        return

    mark(job_id, "processing") # update status
    log.info(f"[{job_id}] processing {input_path}")

    try:
        run_analysis(job_id, input_path, intermediate_dir, result_path, log_dir)
        fasta_paths = sorted((workspace / "results" / "fastas").glob("final_*.fasta"))
        operon_paths = [
            path for path in (operon_plot_path, operon_metadata_path)
            if path.exists()
        ]
        result_paths = [result_path, final_result_path, *fasta_paths, *operon_paths]
        missing = [path for path in result_paths if not path.exists()]
        if missing:
            raise FileNotFoundError(
                "Expected result file(s) were not created: "
                + ", ".join(str(path) for path in missing)
            )
        archive_path = archive_database_job(job, input_path, result_paths, workspace)
        if archive_path:
            log.info("[%s] archived database contribution at %s", job_id, archive_path)
        push_results(job_id, result_paths)
        log.info(f"[{job_id}] complete")
    except subprocess.CalledProcessError as e:
        mark(job_id, "failed", str(e))
        log.exception("[%s] tool failed: %s", job_id, e)
    except Exception as e:
        mark(job_id, "failed", str(e))
        log.exception("[%s] unexpected error: %s", job_id, e)


if __name__ == '__main__':
    log.info(f"Runner started")
    try:
        jobs = poll()
        log.info(f"Polled — {len(jobs)} job(s) ready")
        for job in jobs:
            process(job)
    except Exception as e:
        log.exception("Poll error: %s", e)
