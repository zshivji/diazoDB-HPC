import logging
from sqlmodel import Session

from app.core.db import PUBLIC_USER_ID, engine, init_db
from app.models import User

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_sentinel_user(session: Session) -> None:
    existing = session.get(User, PUBLIC_USER_ID)
    if existing:
        return
    sentinel = User(
        id=PUBLIC_USER_ID,
        email="public@system.internal",
        hashed_password="",       # can never log in
        is_active=False,          # locked out
        is_superuser=False,
        full_name="Public Submissions",
    )
    session.add(sentinel)
    session.commit()


def init() -> None:
    with Session(engine) as session:
        init_db(session)
        create_sentinel_user(session=session)


def main() -> None:
    logger.info("Creating initial data")
    init()
    logger.info("Initial data created")


if __name__ == "__main__":
    main()
