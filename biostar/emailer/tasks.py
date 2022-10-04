import logging
import os
import string

from django.conf import settings

from biostar.emailer import sender

logger = logging.getLogger("engine")


def clean_name(name):
    """
    Strip special chars from the ``name`` portion of a given mail.
    """
    try:
        # Remove punctuation from name
        table = str.maketrans('', '', string.punctuation)
        name = name.translate(table)
    except Exception as exc:
        logger.error(f"Error cleaning name: {name}, {exc}")

    return name


def send_all():
    """
    Needed when the email backend is django-mailer to send queued emails.
    """

    # No email sending during data migration.
    if settings.DATA_MIGRATION:
        return

    # Queued email exists only when the backend is the django-mailer.
    if settings.EMAIL_BACKEND == "mailer.backend.DbBackend":
        try:
            from mailer import engine
            logger.info(f"sending queued emails")
            engine.send_all()
        except Exception as exc:
            logger.error(f"send_all() error: {exc}")


def send_email(template_name, recipient_list, extra_context={}, name="", from_email=None, subject="Subject",
               mass=False):
    """
    Sends an email using a template.
    """

    if not settings.SEND_MAIL:
        return

    # The sender pattern email.
    patt = settings.FROM_EMAIL_PATTERN

    # Final sender email
    from_email = from_email or settings.DEFAULT_FROM_EMAIL
    name = clean_name(name)
    from_email = patt % (name, from_email)

    # Test the templates exists
    if os.path.isfile(template_name):
        logger.error(f"Missing template: {template_name}")
        return False
    try:
        # Generate emails.
        logger.info(f"sending email from={from_email} recipient_list={recipient_list} template={template_name}")

        # The email template instance
        email = sender.EmailTemplate(template_name)

        # Default context added to each template.
        port = f":{settings.HTTP_PORT}"if settings.HTTP_PORT else ""

        context = dict(domain=settings.SITE_DOMAIN, protocol=settings.PROTOCOL,
                       port=port, name=settings.SITE_NAME, subject=subject)

        # Additional context added to the template.
        context.update(extra_context)

        # Generate and send the email.
        if mass:
            email.send_mass(context=context, from_email=from_email, recipient_list=recipient_list)
        else:
            email.send(context=context, from_email=from_email, recipient_list=recipient_list)

        logging.info(f"email sent to recipient_list={recipient_list} ")

    except Exception as exc:
        logger.error(f"send_email error: {exc}")
        return False

    return True
