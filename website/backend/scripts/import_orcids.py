"""Import ORCIDs from the legacy text registry into the contributor table."""

import argparse
import re
from pathlib import Path

from sqlmodel import Session, select

from app.core.db import engine
from app.models import Contributor

ORCID_PATTERN = re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "path",
        nargs="?",
        default="/data/uploads/orcid_ids.txt",
        help="Legacy ORCID text file path",
    )
    args = parser.parse_args()
    path = Path(args.path)
    if not path.exists():
        print(f"No legacy ORCID file found at {path}; nothing to import.")
        return

    orcids = {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if ORCID_PATTERN.fullmatch(line.strip())
    }
    with Session(engine) as session:
        existing = set(session.exec(select(Contributor.orcid)).all())
        added = 0
        for orcid in sorted(orcids - existing):
            session.add(Contributor(orcid=orcid))
            added += 1
        session.commit()

    print(f"Imported {added} ORCID(s); skipped {len(orcids) - added} existing record(s).")


if __name__ == "__main__":
    main()
