"""
Tests for email service — SMTP is mocked.
Run: pytest tests/services/test_email.py -v
"""
from unittest.mock import MagicMock, patch

from app.core.config import settings


def test_send_result_email_calls_smtp(monkeypatch):
    monkeypatch.setattr(settings, "SMTP_HOST", "smtp.test.com")
    monkeypatch.setattr(settings, "SMTP_USER", "user@test.com")
    monkeypatch.setattr(settings, "SMTP_PASSWORD", "secret")
    monkeypatch.setattr(settings, "EMAILS_FROM_EMAIL", "noreply@lab.edu")

    with patch("smtplib.SMTP") as mock_smtp_cls:
        mock_smtp = MagicMock()
        mock_smtp_cls.return_value.__enter__ = lambda s: mock_smtp
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        from app.services.email import send_result_email
        send_result_email(
            to="user@university.edu",
            job_id="abc-123",
            filename="output.csv",
            content_type="text/csv",
            data=b"seq,score\nATCG,0.99",
        )

    mock_smtp.sendmail.assert_called_once()
    args = mock_smtp.sendmail.call_args[0]
    assert args[1] == "user@university.edu"   # recipient
