"""use bigint for job upload sizes

Revision ID: a4f2d8c9b6e3
Revises: 9d3e7f1a2b4c
"""

from alembic import op
import sqlalchemy as sa


revision = "a4f2d8c9b6e3"
down_revision = "9d3e7f1a2b4c"
branch_labels = None
depends_on = None


def upgrade():
    op.alter_column(
        "job",
        "file_size_bytes",
        existing_type=sa.Integer(),
        type_=sa.BigInteger(),
        existing_nullable=True,
    )
    op.alter_column(
        "job",
        "bytes_received",
        existing_type=sa.Integer(),
        type_=sa.BigInteger(),
        existing_nullable=False,
    )


def downgrade():
    op.alter_column(
        "job",
        "bytes_received",
        existing_type=sa.BigInteger(),
        type_=sa.Integer(),
        existing_nullable=False,
    )
    op.alter_column(
        "job",
        "file_size_bytes",
        existing_type=sa.BigInteger(),
        type_=sa.Integer(),
        existing_nullable=True,
    )
