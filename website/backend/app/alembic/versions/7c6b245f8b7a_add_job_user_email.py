"""add job user email

Revision ID: 7c6b245f8b7a
Revises: e9b24197f9d2
Create Date: 2026-05-18 00:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
import sqlmodel.sql.sqltypes


# revision identifiers, used by Alembic.
revision = "7c6b245f8b7a"
down_revision = "e9b24197f9d2"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "job",
        sa.Column("user_email", sqlmodel.sql.sqltypes.AutoString(length=255), nullable=True),
    )
    op.create_index(op.f("ix_job_user_email"), "job", ["user_email"], unique=False)


def downgrade():
    op.drop_index(op.f("ix_job_user_email"), table_name="job")
    op.drop_column("job", "user_email")
