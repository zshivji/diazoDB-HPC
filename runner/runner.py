import os, time, base64, subprocess, requests, logging
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

API       = os.environ["API_URL"].rstrip("/")
HEADERS   = {"x-runner-token": os.environ["RUNNER_SECRET"]}
JOB_BASE  = Path(os.environ["HPC_JOB_BASE"])
INTERVAL  = int(os.environ.get("POLL_INTERVAL_SECONDS", 30))


def poll() -> list[dict]:
    r = requests.get(f"{API}/api/v1/runner/jobs", headers=HEADERS, timeout=15)
    r.raise_for_status()
    return r.json()


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


def push_result(job_id: str, result_path: Path) -> None:
    content_type = (
        "text/csv" if result_path.suffix == ".csv" else "text/html"
    )
    data = base64.b64encode(result_path.read_bytes()).decode()
    r = requests.post(
        f"{API}/api/v1/runner/jobs/{job_id}/result",
        json={"filename": result_path.name, "content_type": content_type, "data_base64": data},
        headers=HEADERS,
        timeout=60,
    )
    r.raise_for_status()


def process(job: dict) -> None:
    job_id   = job["id"]
    hpc_path = Path(job["hpc_path"])
    out_dir  = hpc_path.parent
    out_file = out_dir / "output.csv"

    if not hpc_path.exists():
        log.error(f"[{job_id}] input not found: {hpc_path}")
        mark(job_id, "failed", "Input file not found on HPC")
        return

    mark(job_id, "processing")
    log.info(f"[{job_id}] processing {hpc_path}")

    try:
        # ── Replace this command with your actual analysis tool ──
        subprocess.run(
            ["your_tool", str(hpc_path), "--output", str(out_file)],
            check=True,
            timeout=7200,  # 2-hour hard limit
        )
        push_result(job_id, out_file)
        log.info(f"[{job_id}] complete")
    except subprocess.CalledProcessError as e:
        mark(job_id, "failed", str(e))
        log.error(f"[{job_id}] tool failed: {e}")
    except Exception as e:
        mark(job_id, "failed", str(e))
        log.error(f"[{job_id}] unexpected error: {e}")


if __name__ == '__main__':
    log.info(f"Runner started — polling every {INTERVAL}s")
    while True:
        try:
            jobs = poll()
            log.info(f"Polled — {len(jobs)} job(s) ready")
            for job in jobs:
                process(job)
        except Exception as e:
            log.error(f"Poll error: {e}")
        time.sleep(INTERVAL)