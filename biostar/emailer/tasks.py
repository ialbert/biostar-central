import logging
import os

from django.conf import settings

from biostar.emailer import sender

logger = logging.getLogger("biostar")

def send_all():
    """
    Used when the backend is django-mailer.
    """
    if settings.EMAIL_BACKEND == "mailer.backend.DbBackend":
        try:
            from mailer import engine
            logger.info(f"sending queued emails")
            engine.send_all()
        except Exception as exc:
            logger.error(f"send_all() error: {exc}")


def send_email(template_name, email_list, extra_context={}, from_email=None, subject="Subject"):

    # The senders email.
    from_email = from_email or settings.DEFAULT_FROM_EMAIL

    # Test the templates exists
    if os.path.isfile(template_name):
        logger.error(f"Missing template: {template_name}")
        return False
    try:
        # Generate emails.
        logger.info(f"sending email from={from_email} recipient_list={email_list} template={template_name}")

        # The email template instance
        email = sender.EmailTemplate(template_name)

        # Default context added to each template.
        context = dict(domain=settings.SITE_DOMAIN, protocol=settings.PROTOCOL,
                       port=settings.HTTP_PORT, subject=subject, name=settings.SITE_NAME)

        # Additional context added to the template.
        context.update(extra_context)

        # Generate and send the email.
        email.send(context=context, from_email=from_email, recipient_list=email_list)

        logging.info(f"email sent to recipient_list={email_list} ")

        #1/0

    except Exception as exc:
        logger.error(f"send_email error: {exc}")
        return False

    return True
