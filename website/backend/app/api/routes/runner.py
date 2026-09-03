import binascii
import logging
import uuid
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel
from sqlmodel import Session

from app.api.deps import get_db, verify_runner_token
from app.core.config import settings
from app.crud import get_job, get_jobs_for_runner, update_job
from app.models import Job, JobRunnerView, JobStatus, User
from app.services.email import send_failure_email, send_result_email
from app.services.results import save_result_file

log = logging.getLogger(__name__)

router = APIRouter(
    prefix="/runner",
    tags=["runner"],
    dependencies=[Depends(verify_runner_token)],
)


def _job_url(job_id: uuid.UUID) -> str:
    return f"{settings.FRONTEND_HOST.rstrip('/')}/classify?job={job_id}"


def _job_recipient(*, session: Session, job: Job) -> str | None:
    if job.user_email:
        return job.user_email
    owner = session.get(User, job.owner_id)
    return owner.email if owner else None


@router.get("/jobs", response_model=list[JobRunnerView])
async def poll_jobs(session: Session = Depends(get_db)) -> list[dict]:
    """
    Runner calls this every N seconds.
    Returns all unseen ready jobs.
    """
    # Return ready jobs the runner hasn't seen yet, then mark them seen
    ready_jobs = get_jobs_for_runner(session=session)
    result = []
    for job in ready_jobs:
        hpc_path = str(Path(settings.UPLOAD_DIR) / str(job.id) / job.filename)
        result.append({
            "id": job.id,
            "filename": job.filename,
            "hpc_path": hpc_path,
            "use_prodigal": job.use_prodigal,
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

    email_sent = False
    email_status = "not_applicable"
    if status_val == JobStatus.failed:
        recipient = _job_recipient(session=session, job=job)
        email_status = "recipient_missing"
        if recipient and settings.emails_enabled:
            try:
                send_failure_email(
                    to=recipient,
                    job_id=str(job.id),
                    error_message=payload.get("error_message"),
                )
                email_sent = True
                email_status = "sent"
                log.info("Sent failure email for job %s", job.id)
            except Exception as exc:
                email_status = f"failed: {type(exc).__name__}"
                log.exception("Failed to send failure email for job %s", job.id)
        elif recipient:
            email_status = "smtp_not_configured"
            log.warning(
                "Email is not configured; skipping failure email for job %s",
                job.id,
            )
        else:
            log.warning("No failure-email recipient is available for job %s", job.id)

        update_job(session=session, job=job, email_status=email_status)

    return {"ok": True, "email_sent": email_sent, "email_status": email_status}


@router.get("/jobs/{job_id}/input")
def download_job_input(
    job_id: uuid.UUID,
    session: Session = Depends(get_db),
) -> FileResponse:
    """
    Runner downloads submitted input directly from the API.
    This supports HPC execution without inbound SSH/SCP to the cluster.
    """
    job = get_job(session=session, job_id=job_id)
    if not job:
        raise HTTPException(status_code=404)

    input_path = Path(settings.UPLOAD_DIR) / str(job.id) / job.filename
    if not input_path.exists() or not input_path.is_file():
        raise HTTPException(status_code=404, detail="Input file not found")

    return FileResponse(
        path=input_path,
        filename=job.filename,
        media_type="application/octet-stream",
    )


class ResultPayload(BaseModel):
    # The singular fields preserve compatibility with older runners.
    filename: str | None = None
    content_type: str | None = None
    data_base64: str | None = None
    files: list[dict] | None = None


@router.post("/jobs/{job_id}/result")
def post_result(
    job_id: uuid.UUID,
    payload: ResultPayload,
    session: Session = Depends(get_db),
) -> dict:
    """
    Runner posts the result file here once processing is done.
    API stores the result, emails the submitter when email is configured,
    and marks the job complete.
    """
    import base64

    job = get_job(session=session, job_id=job_id)
    if not job:
        raise HTTPException(status_code=404)

    if payload.files is not None:
        file_payloads = payload.files
    elif payload.filename and payload.content_type and payload.data_base64:
        file_payloads = [{
            "filename": payload.filename,
            "content_type": payload.content_type,
            "data_base64": payload.data_base64,
        }]
    else:
        raise HTTPException(status_code=400, detail="No result files provided")
    if not file_payloads:
        raise HTTPException(status_code=400, detail="No result files provided")

    attachments: list[tuple[str, str, bytes]] = []
    try:
        for file_payload in file_payloads:
            filename = file_payload["filename"]
            content_type = file_payload["content_type"]
            data = base64.b64decode(file_payload["data_base64"], validate=True)
            save_result_file(job_id=job.id, filename=filename, data=data)
            attachments.append((filename, content_type, data))
    except (KeyError, TypeError, binascii.Error, ValueError):
        raise HTTPException(status_code=400, detail="Invalid result payload")

    recipient = _job_recipient(session=session, job=job)

    email_sent = False
    email_status = "recipient_missing"
    if recipient and settings.emails_enabled:
        try:
            send_result_email(
                to=recipient,
                job_id=str(job.id),
                filename=attachments[0][0],
                content_type=attachments[0][1],
                data=attachments[0][2],
                attachments=attachments,
                job_url=_job_url(job.id),
            )
            email_sent = True
            email_status = "sent"
            log.info("Sent result email for job %s", job.id)
        except Exception as exc:
            email_status = f"failed: {type(exc).__name__}"
            log.exception("Failed to send result email for job %s", job.id)
    elif recipient:
        email_status = "smtp_not_configured"
        log.warning("Email is not configured; skipping result email for job %s", job.id)
    else:
        log.warning("No result-email recipient is available for job %s", job.id)

    admin_recipient = settings.JOB_NOTIFICATION_EMAIL
    if (
        admin_recipient
        and str(admin_recipient).lower() != str(recipient).lower()
        and settings.emails_enabled
    ):
        try:
            send_result_email(
                to=admin_recipient,
                job_id=str(job.id),
                filename=attachments[0][0],
                content_type=attachments[0][1],
                data=attachments[0][2],
                attachments=attachments,
                job_url=_job_url(job.id),
            )
            log.info("Sent completion notification for job %s", job.id)
        except Exception:
            log.exception("Failed to send completion notification for job %s", job.id)

    update_job(
        session=session,
        job=job,
        status=JobStatus.complete,
        result_filename=next(
            (filename for filename, _, _ in attachments if filename == "nif_clusters.csv"),
            attachments[0][0],
        ),
        email_status=email_status,
    )
    return {
        "ok": True,
        "email_sent": email_sent,
        "email_status": email_status,
    }
