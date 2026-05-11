.PHONY: reset migrate

reset:
	docker compose down -v
	rm -f backend/app/alembic/versions/*.py
	docker compose up db -d
	sleep 5
	cd backend && uv run alembic revision --autogenerate -m "initial"
	docker compose build backend prestart
	docker compose up -d

migrate:
	docker compose exec backend bash scripts/generate-client.sh
	docker compose build frontend backend prestart
	docker compose up -d
