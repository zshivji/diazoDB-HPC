"""
backend/app/api/routes/jobs_public.py

Public endpoints for the classify page — no JWT auth required.
Separate from /api/v1/jobs/ which requires a logged-in user.
"""
import logging
import uuid
from pathlib import Path

import aiofiles
from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.responses import FileResponse
from pydantic import EmailStr
from sqlmodel import Session, SQLModel

from app.api.deps import get_db
from app.core.config import settings
from app.crud import get_job, update_job
from app.models import Job, JobStatus
from app.services.results import get_result_path, safe_result_filename

log = logging.getLogger(__name__)
router = APIRouter(prefix="/classify", tags=["classify"])


class PublicJobCreate(SQLModel):
    user_email: EmailStr
    filename: str
    file_size_bytes: int
    use_prodigal: bool = False
    sequences: str | None = None  # File content sent directly in JSON


class PublicJobPublic(SQLModel):
    id: uuid.UUID
    filename: str
    file_size_bytes: int | None
    status: JobStatus
    result_filename: str | None = None

def _parse_content_range(content_range: str) -> tuple[int, int, int]:
    try:
        parts = content_range.replace("bytes ", "").split("/")
        start, end = [int(x) for x in parts[0].split("-")]
        total = int(parts[1])
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid Content-Range header")

    return start, end, total


@router.post("/", response_model=PublicJobPublic)
async def create_public_job(
    job_in: PublicJobCreate,
    session: Session = Depends(get_db),
) -> Job:
    """
    Public submission endpoint for the classify page.
    No auth required — user identified by email only.
    Creates a job record and upload directory.
    If file content is provided in JSON, writes it and triggers transfer immediately.
    """
    inline_sequences = job_in.sequences
    if inline_sequences is not None and not inline_sequences.strip():
        raise HTTPException(status_code=400, detail="No sequences provided")

    # Use a sentinel UUID for public jobs (no real owner)
    PUBLIC_OWNER_ID = uuid.UUID("00000000-0000-0000-0000-000000000000")

    job = Job(
        owner_id=PUBLIC_OWNER_ID,
        filename=job_in.filename,
        file_size_bytes=job_in.file_size_bytes,
        status=JobStatus.ready if inline_sequences is not None else JobStatus.created,
        user_email=job_in.user_email,
        use_prodigal=job_in.use_prodigal,
    )
    session.add(job)
    session.commit()
    session.refresh(job)

    # Create upload directory
    job_dir = Path(settings.UPLOAD_DIR) / str(job.id)
    job_dir.mkdir(parents=True, exist_ok=True)

    # If file content is included in the create request, write it now.
    # Otherwise the caller can upload chunks to PATCH /classify/{job_id}/upload.
    if inline_sequences is not None:
        file_path = job_dir / job_in.filename
        file_path.write_text(inline_sequences, encoding="utf-8")
        bytes_written = len(inline_sequences.encode("utf-8"))
        log.info(f"[PUBLIC JOB] {job.id} Wrote {bytes_written} bytes to {file_path}")

        update_job(
            session=session,
            job=job,
            bytes_received=bytes_written,
            status=JobStatus.ready,
        )

    return job


@router.get("/{job_id}", response_model=PublicJobPublic)
def get_public_job_status(
    job_id: uuid.UUID,
    session: Session = Depends(get_db),
) -> Job:
    """
    Public status polling endpoint.
    Returns just enough info for the waiting page to know when to advance.
    """
    job = get_job(session=session, job_id=job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


@router.get("/{job_id}/results/{filename}")
def download_public_job_result(
    job_id: uuid.UUID,
    filename: str,
    session: Session = Depends(get_db),
) -> FileResponse:
    """
    Download a completed result file by job ID.
    Intended for links sent in completion emails.
    """
    job = get_job(session=session, job_id=job_id)
    if not job or job.status != JobStatus.complete:
        raise HTTPException(status_code=404, detail="Result not found")

    try:
        result_path = get_result_path(job_id, filename)
        safe_name = safe_result_filename(filename)
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid result filename")

    if not result_path.exists() or not result_path.is_file():
        raise HTTPException(status_code=404, detail="Result not found")

    return FileResponse(
        path=result_path,
        filename=safe_name,
        media_type="application/octet-stream",
    )


@router.patch("/{job_id}/upload", status_code=200)
async def upload_public_file_chunk(
    job_id: uuid.UUID,
    request: Request,
    session: Session = Depends(get_db),
) -> dict:
    """
    Public file upload endpoint (no auth required).
    Receives chunks with Content-Range header, writes to disk.
    On completion, marks the job ready for the runner.
    """
    job = get_job(session=session, job_id=job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    content_range = request.headers.get("content-range", "")
    log.info(f"[PUBLIC UPLOAD] {job_id} Content-Range: {content_range}")
    start, end, total = _parse_content_range(content_range)

    dest = Path(settings.UPLOAD_DIR) / str(job.id) / job.filename
    dest.parent.mkdir(parents=True, exist_ok=True)

    chunk = await request.body()
    async with aiofiles.open(dest, "r+b" if dest.exists() else "wb") as f:
        await f.seek(start)
        await f.write(chunk)

    new_received = end + 1
    is_complete = new_received >= total

    log.info(f"[PUBLIC UPLOAD] {job_id} bytes_received={new_received}/{total}, complete={is_complete}")

    update_job(
        session=session,
        job=job,
        bytes_received=new_received,
        status=JobStatus.ready if is_complete else JobStatus.uploading,
        file_size_bytes=total,
    )

    return {"bytes_received": new_received, "complete": is_complete}
