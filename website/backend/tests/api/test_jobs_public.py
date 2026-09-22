"""Tests for public classify submissions."""

from fastapi.testclient import TestClient
from sqlmodel import Session, select

from app.models import Contributor, Job


def test_public_job_records_database_consent_and_orcid(
    client: TestClient, db: Session
):
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
    assert db.exec(
        select(Contributor).where(Contributor.orcid == "0000-0000-0000-0001")
    ).first() is not None


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


def test_public_contributors_returns_unique_orcids(
    client: TestClient, db: Session
):
    db.add_all(
        [
            Contributor(
                orcid="0000-0000-0000-0002",
            ),
            Contributor(
                orcid="0000-0000-0000-0001",
            ),
        ]
    )
    db.commit()

    response = client.get("/api/v1/classify/contributors")

    assert response.status_code == 200
    assert response.json() == ["0000-0000-0000-0001", "0000-0000-0000-0002"]