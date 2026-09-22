"""Tests for public classify submissions."""

from fastapi.testclient import TestClient
from sqlmodel import Session

from app.models import Job


def test_public_job_records_database_consent_and_orcid(
    client: TestClient, db: Session, monkeypatch, tmp_path
):
    orcid_file = tmp_path / "orcid_ids.txt"
    monkeypatch.setattr("app.core.config.settings.ORCID_IDS_FILE", str(orcid_file))

    response = client.post(
        "/api/v1/classify/",
        json={
            "user_email": "submitter@example.edu",
            "filename": "sample.fasta",
            "file_size_bytes": 12,
            "include_in_database": True,
            "orcid": "0000-0000-0000-0001",
            "sequences": ">seq1\nMKK\n",
        },
    )

    assert response.status_code == 200
    job = db.get(Job, response.json()["id"])
    assert job is not None
    assert job.include_in_database is True
    assert job.orcid == "0000-0000-0000-0001"
    assert orcid_file.read_text(encoding="utf-8") == "0000-0000-0000-0001\n"


def test_public_job_rejects_invalid_orcid(client: TestClient):
    response = client.post(
        "/api/v1/classify/",
        json={
            "user_email": "submitter@example.edu",
            "filename": "sample.fasta",
            "file_size_bytes": 12,
            "orcid": "not-an-orcid",
            "sequences": ">seq1\nMKK\n",
        },
    )

    assert response.status_code == 422


def test_public_contributors_returns_unique_orcids(client: TestClient, monkeypatch, tmp_path):
    orcid_file = tmp_path / "orcid_ids.txt"
    orcid_file.write_text(
        "0000-0000-0000-0002\ninvalid\n0000-0000-0000-0001\n0000-0000-0000-0002\n",
        encoding="utf-8",
    )
    monkeypatch.setattr("app.core.config.settings.ORCID_IDS_FILE", str(orcid_file))

    response = client.get("/api/v1/classify/contributors")

    assert response.status_code == 200
    assert response.json() == ["0000-0000-0000-0001", "0000-0000-0000-0002"]