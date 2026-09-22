"""Public database records assembled from contributor submissions."""

from fastapi import APIRouter, Depends
from sqlmodel import Session, select

from app.api.deps import get_db
from app.models import UserDatabaseRecord

router = APIRouter(prefix="/database", tags=["database"])


@router.get("/user-records", response_model=list[dict[str, str]])
def get_user_records(session: Session = Depends(get_db)) -> list[dict[str, str]]:
    records = session.exec(
        select(UserDatabaseRecord).order_by(
            UserDatabaseRecord.created_at, UserDatabaseRecord.row_number
        )
    ).all()
    return [{**record.data, "Provenance": record.provenance} for record in records]