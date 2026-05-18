"""
Thin wrapper around the Globus Transfer REST API.
Uses client-credentials grant — no user login needed.

Set GLOBUS_MOCK=true in .env to skip real Globus calls during local dev.
Set COPY_TO_LOCAL=true to copy uploaded files to LOCAL_COPY_DIR instead of using Globus.
"""
import logging
import shutil
from pathlib import Path

import httpx

from app.core.config import settings

log = logging.getLogger(__name__)

_AUTH_URL = "https://auth.globus.org/v2/oauth2/token"
_TRANSFER_URL = "https://transfer.api.globus.org/v0.10"


async def _get_token() -> str:
    async with httpx.AsyncClient() as client:
        r = await client.post(
            _AUTH_URL,
            data={
                "grant_type": "client_credentials",
                "scope": "urn:globus:auth:scope:transfer.api.globus.org:all",
            },
            auth=(settings.GLOBUS_CLIENT_ID, settings.GLOBUS_CLIENT_SECRET),
        )
        r.raise_for_status()
        return r.json()["access_token"]


async def start_transfer(job_id: str) -> str:
    """
    Submits a Globus transfer task for job_id.
    Returns the Globus task_id (store this to poll status later).

    If COPY_TO_LOCAL=true, copies uploaded files to LOCAL_COPY_DIR instead.
    """
    if settings.COPY_TO_LOCAL:
        src = Path(settings.UPLOAD_DIR) / job_id
        dest = Path(settings.LOCAL_COPY_DIR) / job_id
        log.info("[COPY_TO_LOCAL] Attempting to copy %s to %s", src, dest)

        if not src.exists():
            raise FileNotFoundError(f"Upload source not found: {src}")

        src_resolved = src.resolve()
        dest_resolved = dest.resolve() if dest.exists() else dest
        if src_resolved == dest_resolved:
            log.info("[COPY_TO_LOCAL] Source and destination are the same path")
            return f"local-task-{job_id}"

        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(src, dest)
        log.info("[COPY_TO_LOCAL] Copied %s to %s", src, dest)
        return f"local-task-{job_id}"

    if settings.GLOBUS_MOCK:
        log.info(f"[GLOBUS MOCK] Would transfer job {job_id} to HPC")
        return f"mock-task-{job_id}"


    token = await _get_token()
    headers = {"Authorization": f"Bearer {token}"}

    async with httpx.AsyncClient() as client:
        # Activate endpoints — no-op if already active
        for ep in [settings.GLOBUS_SOURCE_ENDPOINT, settings.GLOBUS_DEST_ENDPOINT]:
            await client.post(
                f"{_TRANSFER_URL}/endpoint/{ep}/autoactivate",
                headers=headers,
            )

        payload = {
            "DATA_TYPE": "transfer",
            "source_endpoint": settings.GLOBUS_SOURCE_ENDPOINT,
            "destination_endpoint": settings.GLOBUS_DEST_ENDPOINT,
            "DATA": [{
                "DATA_TYPE": "transfer_item",
                "source_path": f"{settings.UPLOAD_DIR}/{job_id}/",
                "destination_path": f"{settings.GLOBUS_DEST_BASE_PATH}/{job_id}/",
                "recursive": True,
            }],
            "notify_on_succeeded": False,
            "store_base_path_info": False,
        }

        r = await client.post(
            f"{_TRANSFER_URL}/transfer",
            json=payload,
            headers=headers,
        )
        r.raise_for_status()
        return r.json()["task_id"]


async def get_transfer_status(task_id: str) -> str:
    """
    Returns one of: ACTIVE | SUCCEEDED | FAILED | INACTIVE
    Mock/local tasks always return SUCCEEDED
    so the runner poll flow advances normally in local dev.
    """
    if settings.COPY_TO_LOCAL or settings.GLOBUS_MOCK or task_id.startswith("mock-task-") or task_id.startswith("local-task-"):
        return "SUCCEEDED"

    token = await _get_token()
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_TRANSFER_URL}/task/{task_id}",
            headers={"Authorization": f"Bearer {token}"},
        )
        r.raise_for_status()
        return r.json()["status"]
