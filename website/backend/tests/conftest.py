import uuid
from collections.abc import Generator

import pytest
from fastapi.testclient import TestClient
from sqlmodel import Session, delete, select

from app.core.config import settings
from app.core.db import engine, init_db
from app.crud import create_job
from app.main import app
from app.models import Item, Job, JobCreate, User
from tests.utils.user import authentication_token_from_email
from tests.utils.utils import get_superuser_token_headers


def _clear_test_data(session: Session) -> None:
    session.execute(delete(Job))
    session.execute(delete(Item))
    session.execute(
        delete(User).where(
            User.id != uuid.UUID("00000000-0000-0000-0000-000000000000")
        )
    )
    session.commit()


@pytest.fixture(scope="session", autouse=True)
def db() -> Generator[Session, None, None]:
    with Session(engine) as session:
        _clear_test_data(session)
        init_db(session)
        yield session
        _clear_test_data(session)


@pytest.fixture(scope="module")
def client() -> Generator[TestClient, None, None]:
    with TestClient(app) as c:
        yield c


@pytest.fixture(scope="module")
def superuser_token_headers(client: TestClient) -> dict[str, str]:
    return get_superuser_token_headers(client)


@pytest.fixture(scope="module")
def normal_user_token_headers(client: TestClient, db: Session) -> dict[str, str]:
    return authentication_token_from_email(
        client=client, email=settings.EMAIL_TEST_USER, db=db
    )


def _clear_jobs_and_items(session: Session) -> None:
    session.execute(delete(Job))
    session.execute(delete(Item))
    session.commit()


@pytest.fixture(autouse=True)
def clean_records(db) -> Generator[None, None, None]:
    _clear_jobs_and_items(db)
    yield
    _clear_jobs_and_items(db)


@pytest.fixture(autouse=True)
def upload_dir(tmp_path, monkeypatch) -> None:
    path = tmp_path / "uploads"
    path.mkdir()
    monkeypatch.setattr(settings, "UPLOAD_DIR", str(path))


@pytest.fixture
def job(db, request) -> Job:
    """A freshly created job owned by the normal test user."""
    request.getfixturevalue("normal_user_token_headers")
    user = db.exec(select(User).where(User.email == settings.EMAIL_TEST_USER)).first()
    assert user
    return create_job(
        session=db,
        job_in=JobCreate(filename="test.fasta", file_size_bytes=1024),
        owner_id=user.id,
    )


@pytest.fixture
def runner_headers() -> dict:
    """Auth header the HPC runner sends."""
    return {"x-runner-token": settings.RUNNER_SECRET}
