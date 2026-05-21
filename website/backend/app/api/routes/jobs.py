import logging
import uuid
from pathlib import Path

import aiofiles
from fastapi import APIRouter, Depends, HTTPException, Request
from sqlmodel import Session

from app.api.deps import CurrentUser, get_db
from app.core.config import settings
from app.crud import create_job, get_job, update_job
from app.models import Job, JobCreate, JobPublic, JobStatus
from app.services.globus import start_transfer

router = APIRouter(prefix="/jobs", tags=["jobs"])
log = logging.getLogger(__name__)


def _completed_transfer_status() -> JobStatus:
    return JobStatus.ready if settings.COPY_TO_LOCAL or settings.GLOBUS_MOCK or settings.RUNNER_PULL_ONLY else JobStatus.transferring


def _parse_content_range(content_range: str) -> tuple[int, int, int]:
    try:
        parts = content_range.replace("bytes ", "").split("/")
        start, end = [int(x) for x in parts[0].split("-")]
        total = int(parts[1])
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid Content-Range header")

    return start, end, total


@router.post("/", response_model=JobPublic)
async def create_new_job(
    job_in: JobCreate,
    current_user: CurrentUser,
    session: Session = Depends(get_db),
) -> Job:
    """
    Step 1 — create a job record and an empty upload directory.
    The client receives the job ID and begins uploading chunks.
    """
    job = create_job(session=session, job_in=job_in, owner_id=current_user.id)

    # Create the directory that will hold the uploaded file
    job_dir = Path(settings.UPLOAD_DIR) / str(job.id)
    job_dir.mkdir(parents=True, exist_ok=True)

    return job


@router.patch("/{job_id}/upload", status_code=200)
async def upload_chunk(
    job_id: uuid.UUID,
    request: Request,
    current_user: CurrentUser,
    session: Session = Depends(get_db),
) -> dict:
    """
    Step 2 — receive a file chunk.

    Client sends:
      Content-Range: bytes <start>-<end>/<total>
      Body: raw bytes for this chunk

    Supports resume: if the connection drops, the client re-sends
    from the last confirmed offset (GET /{job_id} returns bytes_received).
    """

    job = get_job(session=session, job_id=job_id)
    if not job or job.owner_id != current_user.id:
        raise HTTPException(status_code=404)

    content_range = request.headers.get("content-range", "")
    log.info(f"[UPLOAD] {job_id} Content-Range: {content_range}")
    start, end, total = _parse_content_range(content_range)

    dest = Path(settings.UPLOAD_DIR) / str(job.id) / job.filename
    dest.parent.mkdir(parents=True, exist_ok=True)

    chunk = await request.body()
    async with aiofiles.open(dest, "r+b" if dest.exists() else "wb") as f:
        await f.seek(start)
        await f.write(chunk)

    new_received = end + 1
    is_complete = new_received >= total

    log.info(f"[UPLOAD] {job_id} bytes_received={new_received}/{total}, complete={is_complete}")

    update_job(
        session=session,
        job=job,
        bytes_received=new_received,
        status=JobStatus.uploading if not is_complete else _completed_transfer_status(),
        file_size_bytes=total,
    )

    if is_complete:
        if settings.RUNNER_PULL_ONLY:
            update_job(
                session=session,
                job=job,
                status=JobStatus.ready,
            )
        else:
            log.info(f"[UPLOAD] {job_id} COMPLETE - calling start_transfer()")
            task_id = await start_transfer(str(job.id))
            log.info(f"[UPLOAD] {job_id} start_transfer() returned task_id={task_id}")

            update_job(
                session=session,
                job=job,
                globus_task_id=task_id,
                status=JobStatus.transferring,
            )

    return {"bytes_received": new_received, "complete": is_complete}


@router.get("/{job_id}", response_model=JobPublic)
def get_job_status(
    job_id: uuid.UUID,
    current_user: CurrentUser,
    session: Session = Depends(get_db),
) -> Job:
    """
    Returns current job status + bytes_received.
    Frontend polls this to show progress and resume uploads.
    """
    job = get_job(session=session, job_id=job_id)
    if not job or job.owner_id != current_user.id:
        raise HTTPException(status_code=404)
    return job
