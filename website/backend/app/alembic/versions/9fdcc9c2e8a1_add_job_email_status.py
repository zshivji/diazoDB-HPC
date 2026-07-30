"""add job email status

Revision ID: 9fdcc9c2e8a1
Revises: 4d8e9ab2c731
Create Date: 2026-07-30 00:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
import sqlmodel.sql.sqltypes


# revision identifiers, used by Alembic.
revision = "9fdcc9c2e8a1"
down_revision = "4d8e9ab2c731"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "job",
        sa.Column("email_status", sqlmodel.sql.sqltypes.AutoString(), nullable=True),
    )


def downgrade():
    op.drop_column("job", "email_status")
