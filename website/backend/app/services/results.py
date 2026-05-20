import uuid
from pathlib import Path
from urllib.parse import quote

from app.core.config import settings


def safe_result_filename(filename: str) -> str:
    safe_name = Path(filename.replace("\\", "/")).name
    if not safe_name or safe_name in {".", ".."}:
        raise ValueError("Invalid result filename")
    return safe_name


def get_result_path(job_id: uuid.UUID, filename: str) -> Path:
    return Path(settings.UPLOAD_DIR) / str(job_id) / "results" / safe_result_filename(filename)


def save_result_file(*, job_id: uuid.UUID, filename: str, data: bytes) -> Path:
    path = get_result_path(job_id, filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return path


def build_result_download_url(*, job_id: uuid.UUID, filename: str) -> str | None:
    base_url = settings.BACKEND_PUBLIC_URL
    if not base_url:
        return None
    safe_name = quote(safe_result_filename(filename))
    return f"{base_url.rstrip('/')}{settings.API_V1_STR}/classify/{job_id}/results/{safe_name}"
