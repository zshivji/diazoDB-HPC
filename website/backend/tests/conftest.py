from collections.abc import Generator
import uuid
import pytest
from fastapi.testclient import TestClient
from sqlmodel import Session, delete, select

from app.core.config import settings
from app.core.db import engine, init_db
from app.main import app
from app.models import Item, Job, JobCreate, JobStatus, User
from app.crud import create_job
from tests.utils.user import authentication_token_from_email
from tests.utils.utils import get_superuser_token_headers


@pytest.fixture(scope="session", autouse=True)
def db() -> Generator[Session, None, None]:
    with Session(engine) as session:
        init_db(session)
        yield session
        # Delete in FK-safe order: jobs first, then items, then users
        session.execute(delete(Job))
        session.execute(delete(Item))
        # Delete all users EXCEPT the sentinel (it has no test data to clean)
        session.execute(
            delete(User).where(
                User.id != uuid.UUID("00000000-0000-0000-0000-000000000000")
            )
        )
        session.commit()


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


@pytest.fixture
def job(db) -> Job:
    """
    A freshly created job owned by the superuser (always exists after init_db).
    Avoids depending on a normal_user fixture that doesn't exist in this template.
    """
    superuser = db.exec(
        select(User).where(User.email == settings.FIRST_SUPERUSER)
    ).first()
    return create_job(
        session=db,
        job_in=JobCreate(filename="test.fasta", file_size_bytes=1024),
        owner_id=superuser.id,
    )


@pytest.fixture
def runner_headers() -> dict:
    """Auth header the HPC runner sends."""
    return {"x-runner-token": settings.RUNNER_SECRET}