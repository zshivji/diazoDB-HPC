import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from app.core.config import settings


def send_result_email(
    *,
    to: str,
    job_id: str,
    filename: str,
    content_type: str,
    data: bytes,
    download_url: str | None = None,
) -> None:
    assert settings.emails_enabled, "no provided configuration for email variables"

    msg = MIMEMultipart()
    msg["From"] = settings.EMAILS_FROM_EMAIL
    msg["To"] = to
    msg["Subject"] = f"{settings.PROJECT_NAME} analysis results are ready"

    body = f"Your job ({job_id}) completed successfully."
    if download_url:
        body += f"\n\nDownload your result here:\n{download_url}"
    body += "\n\nThe result file is also attached to this email."

    msg.attach(MIMEText(body, "plain"))

    part = MIMEBase(*content_type.split("/"))
    part.set_payload(data)
    encoders.encode_base64(part)
    part.add_header("Content-Disposition", f'attachment; filename="{filename}"')
    msg.attach(part)

    smtp_cls = smtplib.SMTP_SSL if settings.SMTP_SSL else smtplib.SMTP
    with smtp_cls(settings.SMTP_HOST, settings.SMTP_PORT) as server:
        if settings.SMTP_TLS and not settings.SMTP_SSL:
            server.starttls()
        if settings.SMTP_USER:
            server.login(settings.SMTP_USER, settings.SMTP_PASSWORD)
        server.sendmail(settings.EMAILS_FROM_EMAIL, to, msg.as_string())
