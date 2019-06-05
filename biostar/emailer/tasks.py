import os
import logging

from biostar.emailer import sender

from django.conf import settings

logger = logging.getLogger("engine")


def send_email(template_name, email_list, extra_context={}, from_email=None, subject="Subject", send=False):
    from mailer.engine import send_all

    from_email = from_email or settings.ADMINS[0][1]

    # Test the templates exists
    if os.path.isfile(template_name):
        logger.error(f"Missing template: {template_name}")
        return False
    try:
        # Generate emails.
        logger.debug(f"Emails from={from_email} to email_list={email_list} using template:{template_name}")

        # The object that parsers the template.
        email = sender.EmailTemplate(template_name)

        # Accumulate the emails into the database.
        context = dict(domain=settings.SITE_DOMAIN, protocol=settings.PROTOCOL, port=settings.HTTP_PORT,
                       subject=subject, name=settings.SITE_NAME)
        context.update(extra_context)
        email.send(context=context, from_email=from_email, recipient_list=email_list)

        if send and email_list:
            # Send the emails.
            send_all()
            logging.info(f"Email has been sent to { len(email_list) } accounts.")

    except Exception as exc:
        logger.error(f"Mailing error : {exc}")
        return False

    return True
