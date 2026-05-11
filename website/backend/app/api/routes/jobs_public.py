"""
backend/app/api/routes/jobs_public.py
 
Public endpoints for the classify page — no JWT auth required.
Separate from /api/v1/jobs/ which requires a logged-in user.
"""
import uuid
from pathlib import Path
from fastapi import APIRouter, Depends, HTTPException
from sqlmodel import Session, SQLModel
 
from app.api.deps import get_db
from app.core.config import settings
from app.models import Job, JobPublic, JobStatus
from app.crud import get_job
 
router = APIRouter(prefix="/classify", tags=["classify"])
 
 
class PublicJobCreate(SQLModel):
    user_email: str
    filename: str
    file_size_bytes: int
    use_prodigal: bool = False
 
 
class PublicJobPublic(SQLModel):
    id: uuid.UUID
    filename: str
    file_size_bytes: int | None
    status: JobStatus
 
 
@router.post("/", response_model=PublicJobPublic)
def create_public_job(
    job_in: PublicJobCreate,
    session: Session = Depends(get_db),
) -> Job:
    """
    Public submission endpoint for the classify page.
    No auth required — user identified by email only.
    Creates a job record and upload directory.
    The frontend then uploads the FASTA content separately.
    """
    # Use a sentinel UUID for public jobs (no real owner)
    PUBLIC_OWNER_ID = uuid.UUID("00000000-0000-0000-0000-000000000000")
 
    job = Job(
        owner_id=PUBLIC_OWNER_ID,
        filename=job_in.filename,
        file_size_bytes=job_in.file_size_bytes,
        status=JobStatus.created,
        # Store user email in error_message field temporarily
        # TODO: add user_email column to Job model in a future migration
        error_message=f"public:{job_in.user_email}",
    )
    session.add(job)
    session.commit()
    session.refresh(job)
 
    # Create upload directory
    job_dir = Path(settings.UPLOAD_DIR) / str(job.id)
    job_dir.mkdir(parents=True, exist_ok=True)
 
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
 