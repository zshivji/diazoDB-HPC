"""
Tests for user-facing job endpoints.
Run: pytest tests/api/test_jobs.py -v
"""
from unittest.mock import AsyncMock, patch

from fastapi.testclient import TestClient

# ── Step 1: job creation ──────────────────────────────────────────────────────

def test_create_job_returns_id(client: TestClient, normal_user_token_headers):
    r = client.post(
        "/api/v1/jobs/",
        json={"filename": "sample.fasta", "file_size_bytes": 1_000_000},
        headers=normal_user_token_headers,
    )
    assert r.status_code == 200
    data = r.json()
    assert "id" in data
    assert data["status"] == "created"
    assert data["filename"] == "sample.fasta"


def test_create_job_requires_auth(client: TestClient):
    r = client.post(
        "/api/v1/jobs/",
        json={"filename": "x.fasta", "file_size_bytes": 100},
    )
    assert r.status_code == 401


# ── Step 2: chunked upload ────────────────────────────────────────────────────

def test_upload_single_chunk(client, normal_user_token_headers, job, tmp_path, monkeypatch):
    """A single chunk that completes the upload should trigger Globus transfer."""
    monkeypatch.setattr("app.core.config.settings.UPLOAD_DIR", str(tmp_path))

    fake_data = b"ATCG" * 256   # 1 KB
    total = len(fake_data)

    with patch("app.api.routes.jobs.start_transfer", new_callable=AsyncMock) as mock_globus:
        mock_globus.return_value = "fake-globus-task-id"

        r = client.patch(
            f"/api/v1/jobs/{job.id}/upload",
            content=fake_data,
            headers={
                **normal_user_token_headers,
                "Content-Range": f"bytes 0-{total - 1}/{total}",
            },
        )

    assert r.status_code == 200
    body = r.json()
    assert body["bytes_received"] == total
    assert body["complete"] is True
    mock_globus.assert_called_once_with(str(job.id))


def test_upload_chunk_bad_content_range(client, normal_user_token_headers, job):
    r = client.patch(
        f"/api/v1/jobs/{job.id}/upload",
        content=b"data",
        headers={**normal_user_token_headers, "Content-Range": "garbage"},
    )
    assert r.status_code == 400


def test_upload_resumes_from_offset(client, normal_user_token_headers, job, tmp_path, monkeypatch):
    """
    Simulate a dropped connection: upload first half, then 'reconnect' and
    GET the job to find bytes_received, then send second half.
    """
    monkeypatch.setattr("app.core.config.settings.UPLOAD_DIR", str(tmp_path))
    (tmp_path / str(job.id)).mkdir(parents=True, exist_ok=True)

    total = 1000
    first_half = b"A" * 500

    with patch("app.api.routes.jobs.start_transfer", new_callable=AsyncMock):
        # First chunk
        r = client.patch(
            f"/api/v1/jobs/{job.id}/upload",
            content=first_half,
            headers={
                **normal_user_token_headers,
                "Content-Range": f"bytes 0-499/{total}",
            },
        )
    assert r.json()["bytes_received"] == 500

    # Client checks progress
    r = client.get(f"/api/v1/jobs/{job.id}", headers=normal_user_token_headers)
    assert r.json()["bytes_received"] == 500

    with patch("app.api.routes.jobs.start_transfer", new_callable=AsyncMock) as m:
        m.return_value = "task-id"
        # Second (final) chunk
        r = client.patch(
            f"/api/v1/jobs/{job.id}/upload",
            content=b"B" * 500,
            headers={
                **normal_user_token_headers,
                "Content-Range": f"bytes 500-999/{total}",
            },
        )
    assert r.json()["complete"] is True


# ── Step 3: job status ────────────────────────────────────────────────────────

def test_get_job_status(client, normal_user_token_headers, job):
    r = client.get(f"/api/v1/jobs/{job.id}", headers=normal_user_token_headers)
    assert r.status_code == 200
    assert r.json()["id"] == str(job.id)


def test_get_job_other_user_is_404(client, superuser_token_headers, job):
    """Superuser cannot see another user's job via the regular endpoint."""
    r = client.get(f"/api/v1/jobs/{job.id}", headers=superuser_token_headers)
    assert r.status_code == 404
