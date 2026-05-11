import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email import encoders
from app.core.config import settings


def send_result_email(
    *,
    to: str,
    job_id: str,
    filename: str,
    content_type: str,
    data: bytes,
) -> None:
    msg = MIMEMultipart()
    msg["From"] = settings.EMAILS_FROM_ADDRESS
    msg["To"] = to
    msg["Subject"] = "Your FASTA analysis is ready"

    msg.attach(MIMEText(
        f"Your job ({job_id}) completed successfully. Results are attached.",
        "plain",
    ))

    part = MIMEBase(*content_type.split("/"))
    part.set_payload(data)
    encoders.encode_base64(part)
    part.add_header("Content-Disposition", f'attachment; filename="{filename}"')
    msg.attach(part)

    with smtplib.SMTP(settings.SMTP_HOST, settings.SMTP_PORT) as server:
        server.starttls()
        server.login(settings.SMTP_USER, settings.SMTP_PASSWORD)
        server.sendmail(settings.EMAILS_FROM_ADDRESS, to, msg.as_string())