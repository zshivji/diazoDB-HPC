import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from app.core.config import settings


def _send_message(msg: MIMEMultipart) -> None:
    smtp_cls = smtplib.SMTP_SSL if settings.SMTP_SSL else smtplib.SMTP
    with smtp_cls(settings.SMTP_HOST, settings.SMTP_PORT) as server:
        if settings.SMTP_TLS and not settings.SMTP_SSL:
            server.starttls()
        if settings.SMTP_USER:
            server.login(settings.SMTP_USER, settings.SMTP_PASSWORD)
        server.sendmail(settings.EMAILS_FROM_EMAIL, msg["To"], msg.as_string())


def send_submission_email(
    *,
    to: str,
    job_id: str,
    submitter_email: str,
    filename: str,
    file_size_bytes: int,
    job_url: str | None = None,
) -> None:
    """Notify the configured testing/admin inbox about a new submission."""
    assert settings.emails_enabled, "no provided configuration for email variables"

    msg = MIMEMultipart()
    msg["From"] = settings.EMAILS_FROM_EMAIL
    msg["To"] = to
    msg["Subject"] = f"{settings.PROJECT_NAME} new job submitted"

    body = (
        f"A new {settings.PROJECT_NAME} job was submitted.\n\n"
        f"Job ID: {job_id}\n"
        f"Submitter email: {submitter_email}\n"
        f"Input file: {filename}\n"
        f"File size: {file_size_bytes} bytes"
    )
    if job_url:
        body += f"\n\nView the job here:\n{job_url}"
    msg.attach(MIMEText(body, "plain"))
    _send_message(msg)


def send_result_email(
    *,
    to: str,
    job_id: str,
    filename: str,
    content_type: str,
    data: bytes,
    download_url: str | None = None,
    job_url: str | None = None,
    attachments: list[tuple[str, str, bytes]] | None = None,
) -> None:
    assert settings.emails_enabled, "no provided configuration for email variables"

    msg = MIMEMultipart()
    msg["From"] = settings.EMAILS_FROM_EMAIL
    msg["To"] = to
    msg["Subject"] = f"{settings.PROJECT_NAME} analysis results are ready"

    body = f"Your job ({job_id}) completed successfully."
    if job_url:
        body += f"\n\nView your job and download the results here:\n{job_url}"
    elif download_url:
        body += f"\n\nDownload your result here:\n{download_url}"

    if attachments is None:
        attachments = [(filename, content_type, data)]
    body += "\n\nThe result files are also attached to this email."

    msg.attach(MIMEText(body, "plain"))

    for attachment_name, attachment_type, attachment_data in attachments:
        part = MIMEBase(*attachment_type.split("/"))
        part.set_payload(attachment_data)
        encoders.encode_base64(part)
        part.add_header(
            "Content-Disposition", f'attachment; filename="{attachment_name}"'
        )
        msg.attach(part)

    _send_message(msg)


def send_failure_email(
    *,
    to: str,
    job_id: str,
    error_message: str | None = None,
) -> None:
    assert settings.emails_enabled, "no provided configuration for email variables"

    msg = MIMEMultipart()
    msg["From"] = settings.EMAILS_FROM_EMAIL
    msg["To"] = to
    msg["Subject"] = f"{settings.PROJECT_NAME} analysis failed"

    body = f"Your job ({job_id}) did not complete successfully."
    if error_message:
        body += f"\n\nError:\n{error_message}"

    msg.attach(MIMEText(body, "plain"))
    _send_message(msg)
