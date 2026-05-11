"""
Tests for Globus service — all HTTP is mocked with pytest-mock.
Run: pytest tests/services/test_globus.py -v
"""
import pytest
import httpx


@pytest.fixture(autouse=True)
def disable_globus_mock(monkeypatch):
    """
    Force GLOBUS_MOCK=False for every test in this file.
    Without this, start_transfer returns early with a mock-task-id
    before any HTTP calls happen, making the patches useless.
    """
    from app.core.config import settings
    monkeypatch.setattr(settings, "GLOBUS_MOCK", False)


@pytest.mark.anyio
async def test_start_transfer_returns_task_id(mocker):
    """start_transfer should POST to Globus and return task_id."""
    token_response = httpx.Response(200, json={"access_token": "fake-token"})
    activate_response = httpx.Response(200, json={"code": "AutoActivated"})
    transfer_response = httpx.Response(200, json={"task_id": "globus-task-abc"})

    responses = iter([
        token_response,
        activate_response,
        activate_response,
        transfer_response,
    ])

    async def mock_send(self, request, **kwargs):
        return next(responses)

    mocker.patch("httpx.AsyncClient.send", mock_send)

    from app.services.globus import start_transfer
    task_id = await start_transfer("test-job-id")
    assert task_id == "globus-task-abc"


@pytest.mark.anyio
async def test_get_transfer_status_returns_succeeded(mocker):
    token_response = httpx.Response(200, json={"access_token": "fake-token"})
    task_response = httpx.Response(200, json={"status": "SUCCEEDED"})

    call_count = {"n": 0}

    async def mock_send(self, request, **kwargs):
        call_count["n"] += 1
        return token_response if call_count["n"] == 1 else task_response

    mocker.patch("httpx.AsyncClient.send", mock_send)

    from app.services.globus import get_transfer_status
    status = await get_transfer_status("some-task-id")
    assert status == "SUCCEEDED"


@pytest.mark.anyio
async def test_globus_raises_on_auth_failure(mocker):
    async def mock_send(self, request, **kwargs):
        return httpx.Response(401, json={"error": "unauthorized"})

    mocker.patch("httpx.AsyncClient.send", mock_send)

    from app.services.globus import start_transfer
    with pytest.raises(httpx.HTTPStatusError):
        await start_transfer("job-id")