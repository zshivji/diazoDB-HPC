"""add user database records

Revision ID: 9d3e7f1a2b4c
Revises: 8c2d9f4a1b6e
"""

from alembic import op
import sqlalchemy as sa


revision = "9d3e7f1a2b4c"
down_revision = "8c2d9f4a1b6e"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "user_database_record",
        sa.Column("id", sa.Uuid(), nullable=False),
        sa.Column("job_id", sa.Uuid(), nullable=False),
        sa.Column("row_number", sa.Integer(), nullable=False),
        sa.Column("provenance", sa.String(length=32), nullable=False),
        sa.Column("orcid", sa.String(length=19), nullable=True),
        sa.Column("data", sa.JSON(), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.ForeignKeyConstraint(["job_id"], ["job.id"]),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("job_id", "row_number"),
    )
    op.create_index(
        "ix_user_database_record_job_id",
        "user_database_record",
        ["job_id"],
        unique=False,
    )
    op.create_index(
        "ix_user_database_record_provenance",
        "user_database_record",
        ["provenance"],
        unique=False,
    )


def downgrade():
    op.drop_index("ix_user_database_record_provenance", table_name="user_database_record")
    op.drop_index("ix_user_database_record_job_id", table_name="user_database_record")
    op.drop_table("user_database_record")