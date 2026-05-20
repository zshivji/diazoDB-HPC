"""add job use prodigal

Revision ID: 4d8e9ab2c731
Revises: 7c6b245f8b7a
Create Date: 2026-05-18 00:00:00.000000

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "4d8e9ab2c731"
down_revision = "7c6b245f8b7a"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "job",
        sa.Column("use_prodigal", sa.Boolean(), nullable=False, server_default=sa.false()),
    )
    op.alter_column("job", "use_prodigal", server_default=None)


def downgrade():
    op.drop_column("job", "use_prodigal")
