#!/bin/sh
# set up nightly backups of diazoDB docker volumes (postgres tables)
set -eu

COMPOSE_DIR=/home/zshivji/diazoDB-HPC/website
BACKUP_DIR=/home/zshivji/diazoDB-HPC/srv-backups/postgres
DATE=$(date +%F_%H%M%S)

mkdir -p "$BACKUP_DIR"

cd "$COMPOSE_DIR"

docker compose exec -T db sh -c \
  'pg_dump -U "$POSTGRES_USER" -d "$POSTGRES_DB" --format=custom' \
  > "$BACKUP_DIR/diazodb-$DATE.dump"

# delete backups after 30 days
find "$BACKUP_DIR" -type f -name 'diazodb-*.dump' -mtime +30 -delete
