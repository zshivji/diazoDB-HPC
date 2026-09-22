"""add database consent and ORCID to jobs

Revision ID: 5b7f1c2d9e31
Revises: 9fdcc9c2e8a1
"""

from alembic import op
import sqlalchemy as sa


revision = "5b7f1c2d9e31"
down_revision = "9fdcc9c2e8a1"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "job",
        sa.Column("include_in_database", sa.Boolean(), nullable=False, server_default=sa.false()),
    )
    op.add_column("job", sa.Column("orcid", sa.String(length=19), nullable=True))
    op.alter_column("job", "include_in_database", server_default=None)


def downgrade():
    op.drop_column("job", "orcid")
    op.drop_column("job", "include_in_database")