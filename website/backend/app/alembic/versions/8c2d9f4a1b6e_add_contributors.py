"""add durable contributor ORCID records

Revision ID: 8c2d9f4a1b6e
Revises: 5b7f1c2d9e31
"""

from alembic import op
import sqlalchemy as sa


revision = "8c2d9f4a1b6e"
down_revision = "5b7f1c2d9e31"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "contributor",
        sa.Column("id", sa.Uuid(), nullable=False),
        sa.Column("orcid", sa.String(length=19), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("orcid"),
    )
    op.create_index(
        "ix_contributor_orcid", "contributor", ["orcid"], unique=False
    )


def downgrade():
    op.drop_index("ix_contributor_orcid", table_name="contributor")
    op.drop_table("contributor")
