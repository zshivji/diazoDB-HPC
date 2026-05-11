from fastapi import APIRouter

from app.api.routes import items, login, private, users, utils, jobs, runner, jobs_public
from app.core.config import settings

api_router = APIRouter()
api_router.include_router(login.router)
api_router.include_router(users.router)
api_router.include_router(utils.router)
api_router.include_router(items.router)
api_router.include_router(jobs.router)
api_router.include_router(runner.router)
api_router.include_router(jobs_public.router)

if settings.ENVIRONMENT == "local":
    api_router.include_router(private.router)