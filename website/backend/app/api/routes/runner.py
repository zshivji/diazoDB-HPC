import uuid

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel
from sqlmodel import Session

from app.api.deps import get_db, verify_runner_token
from app.core.config import settings
from app.crud import get_job, get_jobs_for_runner, update_job
from app.models import Job, JobRunnerView, JobStatus
from app.services.email import send_result_email
from app.services.globus import get_transfer_status

router = APIRouter(
    prefix="/runner",
    tags=["runner"],
    dependencies=[Depends(verify_runner_token)],
)


@router.get("/jobs", response_model=list[JobRunnerView])
async def poll_jobs(session: Session = Depends(get_db)) -> list[dict]:
    """
    Runner calls this every N seconds.
    Lazily checks Globus transfer status for 'transferring' jobs,
    promotes them to 'ready', then returns all unseen ready jobs.
    """
    from sqlmodel import select

    from app.models import JobStatus

    # Check all jobs still mid-transfer
    transferring = session.exec(
        select(Job).where(Job.status == JobStatus.transferring)
    ).all()

    for job in transferring:
        if not job.globus_task_id:
            continue
        globus_status = await get_transfer_status(job.globus_task_id)
        if globus_status == "SUCCEEDED":
            update_job(session=session, job=job, status=JobStatus.ready)
        elif globus_status == "FAILED":
            update_job(session=session, job=job, status=JobStatus.failed,
                       error_message="Globus transfer failed")

    # Return ready jobs the runner hasn't seen yet, then mark them seen
    ready_jobs = get_jobs_for_runner(session=session)
    result = []
    for job in ready_jobs:
        hpc_path = f"{settings.GLOBUS_DEST_BASE_PATH}/{job.id}/{job.filename}"
        result.append({
            "id": job.id,
            "filename": job.filename,
            "hpc_path": hpc_path,
        })
        update_job(session=session, job=job, seen_by_runner=True)

    return result


@router.patch("/jobs/{job_id}")
def update_job_status(
    job_id: uuid.UUID,
    payload: dict,        # {"status": "processing" | "failed"}
    session: Session = Depends(get_db),
) -> dict:
    """Runner calls this to signal it has started or failed a job."""
    job = get_job(session=session, job_id=job_id)
    if not job:
        raise HTTPException(status_code=404)
    status_val = JobStatus(payload["status"])
    kw = {"status": status_val}
    if "error_message" in payload:
        kw["error_message"] = payload["error_message"]
    update_job(session=session, job=job, **kw)
    return {"ok": True}


class ResultPayload(BaseModel):
    filename: str
    content_type: str       # "text/csv" or "text/html"
    data_base64: str        # base64-encoded file content


@router.post("/jobs/{job_id}/result")
def post_result(
    job_id: uuid.UUID,
    payload: ResultPayload,
    session: Session = Depends(get_db),
) -> dict:
    """
    Runner posts the result file here once processing is done.
    API stores result metadata, sends email to job owner, marks complete.
    """
    import base64

    from app.models import User

    job = get_job(session=session, job_id=job_id)
    if not job:
        raise HTTPException(status_code=404)

    data = base64.b64decode(payload.data_base64)

    # Look up owner email
    owner = session.get(User, job.owner_id)
    if owner:
        send_result_email(
            to=owner.email,
            job_id=str(job.id),
            filename=payload.filename,
            content_type=payload.content_type,
            data=data,
        )

    update_job(
        session=session,
        job=job,
        status=JobStatus.complete,
        result_filename=payload.filename,
    )
    return {"ok": True}
