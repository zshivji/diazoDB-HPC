import uuid
import enum
from datetime import datetime, timezone
from typing import Optional

from pydantic import EmailStr
from sqlmodel import Field, Relationship, SQLModel
from sqlalchemy import DateTime


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_datetime_utc() -> datetime:
    return datetime.now(timezone.utc)


# ── Job models ────────────────────────────────────────────────────────────────

class JobStatus(str, enum.Enum):
    created      = "created"
    uploading    = "uploading"
    transferring = "transferring"   # Globus transfer in flight
    ready        = "ready"          # file landed on HPC, runner can pick up
    processing   = "processing"     # runner is working
    complete     = "complete"
    failed       = "failed"


class Job(SQLModel, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    owner_id: uuid.UUID = Field(foreign_key="user.id", nullable=False)
    filename: str
    file_size_bytes: Optional[int] = None
    bytes_received: int = Field(default=0)
    status: JobStatus = Field(default=JobStatus.created)
    globus_task_id: Optional[str] = None
    seen_by_runner: bool = Field(default=False)
    result_filename: Optional[str] = None
    error_message: Optional[str] = None
    created_at: datetime = Field(default_factory=get_datetime_utc)
    updated_at: datetime = Field(default_factory=get_datetime_utc)


class JobCreate(SQLModel):
    filename: str
    file_size_bytes: int


class JobPublic(SQLModel):
    id: uuid.UUID
    filename: str
    file_size_bytes: Optional[int]
    bytes_received: int
    status: JobStatus
    created_at: datetime


class JobRunnerView(SQLModel):
    """Shape the runner receives when polling."""
    id: uuid.UUID
    filename: str
    hpc_path: str   # absolute path on HPC after Globus transfer


# ── User models ───────────────────────────────────────────────────────────────

class UserBase(SQLModel):
    email: EmailStr = Field(unique=True, index=True, max_length=255)
    is_active: bool = True
    is_superuser: bool = False
    full_name: str | None = Field(default=None, max_length=255)


class UserCreate(UserBase):
    password: str = Field(min_length=8, max_length=128)


class UserRegister(SQLModel):
    email: EmailStr = Field(max_length=255)
    password: str = Field(min_length=8, max_length=128)
    full_name: str | None = Field(default=None, max_length=255)


class UserUpdate(UserBase):
    email: EmailStr | None = Field(default=None, max_length=255)  # type: ignore[assignment]
    password: str | None = Field(default=None, min_length=8, max_length=128)


class UserUpdateMe(SQLModel):
    full_name: str | None = Field(default=None, max_length=255)
    email: EmailStr | None = Field(default=None, max_length=255)


class UpdatePassword(SQLModel):
    current_password: str = Field(min_length=8, max_length=128)
    new_password: str = Field(min_length=8, max_length=128)


class User(UserBase, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    hashed_password: str
    created_at: datetime | None = Field(
        default_factory=get_datetime_utc,
        sa_type=DateTime(timezone=True),  # type: ignore
    )
    items: list["Item"] = Relationship(back_populates="owner", cascade_delete=True)


class UserPublic(UserBase):
    id: uuid.UUID
    created_at: datetime | None = None


class UsersPublic(SQLModel):
    data: list[UserPublic]
    count: int


# ── Item models ───────────────────────────────────────────────────────────────

class ItemBase(SQLModel):
    title: str = Field(min_length=1, max_length=255)
    description: str | None = Field(default=None, max_length=255)


class ItemCreate(ItemBase):
    pass


class ItemUpdate(ItemBase):
    title: str | None = Field(default=None, min_length=1, max_length=255)  # type: ignore[assignment]


class Item(ItemBase, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default_factory=get_datetime_utc,
        sa_type=DateTime(timezone=True),  # type: ignore
    )
    owner_id: uuid.UUID = Field(
        foreign_key="user.id", nullable=False, ondelete="CASCADE"
    )
    owner: User | None = Relationship(back_populates="items")


class ItemPublic(ItemBase):
    id: uuid.UUID
    owner_id: uuid.UUID
    created_at: datetime | None = None


class ItemsPublic(SQLModel):
    data: list[ItemPublic]
    count: int


# ── Auth / utility schemas ────────────────────────────────────────────────────

class Message(SQLModel):
    message: str


class Token(SQLModel):
    access_token: str
    token_type: str = "bearer"


class TokenPayload(SQLModel):
    sub: str | None = None


class NewPassword(SQLModel):
    token: str
    new_password: str = Field(min_length=8, max_length=128)