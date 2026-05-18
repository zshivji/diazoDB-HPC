"""
Tests for runner-facing endpoints.
Run: pytest tests/api/test_runner.py -v
"""
import base64
from unittest.mock import AsyncMock, patch

from fastapi.testclient import TestClient

from app.crud import update_job
from app.models import JobStatus

# ── Step 4: runner auth ───────────────────────────────────────────────────────

def test_runner_poll_rejects_missing_token(client: TestClient):
    r = client.get("/api/v1/runner/jobs")
    assert r.status_code == 422   # missing required header


def test_runner_poll_rejects_wrong_token(client: TestClient):
    r = client.get("/api/v1/runner/jobs", headers={"x-runner-token": "wrong"})
    assert r.status_code == 401


# ── Step 5: polling ───────────────────────────────────────────────────────────

def test_runner_poll_empty_when_no_ready_jobs(client, runner_headers):
    r = client.get("/api/v1/runner/jobs", headers=runner_headers)
    assert r.status_code == 200
    assert r.json() == []


def test_runner_poll_returns_ready_jobs(client, runner_headers, job, db):
    """A job in 'ready' state should appear in the poll response."""
    update_job(session=db, job=job, status=JobStatus.ready, globus_task_id=None)

    with patch("app.api.routes.runner.get_transfer_status", new_callable=AsyncMock):
        r = client.get("/api/v1/runner/jobs", headers=runner_headers)

    assert r.status_code == 200
    jobs = r.json()
    assert len(jobs) == 1
    assert jobs[0]["id"] == str(job.id)
    assert "hpc_path" in jobs[0]


def test_runner_poll_promotes_succeeded_globus_transfer(client, runner_headers, job, db):
    """
    A job stuck in 'transferring' should be promoted to 'ready'
    when Globus reports SUCCEEDED.
    """
    update_job(
        session=db, job=job,
        status=JobStatus.transferring,
        globus_task_id="fake-task-id",
    )

    with patch(
        "app.api.routes.runner.get_transfer_status",
        new_callable=AsyncMock,
        return_value="SUCCEEDED",
    ):
        r = client.get("/api/v1/runner/jobs", headers=runner_headers)

    assert r.status_code == 200
    assert len(r.json()) == 1


def test_runner_poll_marks_jobs_seen(client, runner_headers, job, db):
    """After one poll, the same job must not appear in the next poll."""
    update_job(session=db, job=job, status=JobStatus.ready)

    with patch("app.api.routes.runner.get_transfer_status", new_callable=AsyncMock):
        r1 = client.get("/api/v1/runner/jobs", headers=runner_headers)
        r2 = client.get("/api/v1/runner/jobs", headers=runner_headers)

    assert len(r1.json()) == 1
    assert len(r2.json()) == 0   # already seen


# ── Step 6: status update ─────────────────────────────────────────────────────

def test_runner_can_mark_processing(client, runner_headers, job, db):
    update_job(session=db, job=job, status=JobStatus.ready, seen_by_runner=True)
    r = client.patch(
        f"/api/v1/runner/jobs/{job.id}",
        json={"status": "processing"},
        headers=runner_headers,
    )
    assert r.status_code == 200
    assert r.json()["ok"] is True


def test_runner_can_mark_failed(client, runner_headers, job, db):
    update_job(session=db, job=job, status=JobStatus.processing, seen_by_runner=True)
    r = client.patch(
        f"/api/v1/runner/jobs/{job.id}",
        json={"status": "failed", "error_message": "tool crashed"},
        headers=runner_headers,
    )
    assert r.status_code == 200


# ── Step 7: result + email ────────────────────────────────────────────────────

def test_runner_post_result_sends_email(client, runner_headers, job, db):
    update_job(session=db, job=job, status=JobStatus.processing, seen_by_runner=True)

    csv_bytes = b"seq,score\nATCG,0.99\n"
    payload = {
        "filename": "output.csv",
        "content_type": "text/csv",
        "data_base64": base64.b64encode(csv_bytes).decode(),
    }

    with patch("app.api.routes.runner.send_result_email") as mock_email:
        r = client.post(
            f"/api/v1/runner/jobs/{job.id}/result",
            json=payload,
            headers=runner_headers,
        )

    assert r.status_code == 200
    assert r.json()["ok"] is True
    mock_email.assert_called_once()
    call_kwargs = mock_email.call_args.kwargs
    assert call_kwargs["filename"] == "output.csv"
    assert call_kwargs["data"] == csv_bytes
