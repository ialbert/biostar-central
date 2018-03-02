import os
import logging
from django.contrib.sites.models import Site

from biostar.emailer import sender
from mailer.engine import send_all

from django.conf import settings
from .models import  EmailAddress, Subscription

logger = logging.getLogger("engine")


def add_subscription(email, group, name=''):

    # Get the address from the database.
    address = EmailAddress.objects.filter(email=email).first()
    if not address:
        address = EmailAddress.objects.create(name=name, email=email)

    # Fetch the subscriptions if these may exists.
    query = Subscription.objects.filter(group=group, address=address)

    # Drop subscription if it exists.
    query.delete()

    # Create the new subscription.
    Subscription.objects.create(group=group, address=address)


def notify(template_name, email_list, extra_context={},from_email=None, subject="Subject", send=False):

    from_email = from_email or settings.ADMINS[0][1]

    # Test the templates
    if os.path.isfile(template_name):
        logger.error(f"Missing template: {template_name}")
        return False
    try:
        # Generate emails.
        logger.info(f"Emails from={from_email} to email_list={email_list} using template:{template_name}")

        # The object that parsers the template.
        email = sender.EmailTemplate(template_name)

        # This is the context passed to each template.
        site = Site.objects.get_current()

        # Accumulate the emails into the database.
        context= dict(site=site, protocol=settings.PROTOCOL, subject=subject)
        context.update(extra_context)
        email.send(context=context, from_email=from_email, recipient_list=email_list)

        if send:
            # Send the emails.
            send_all()
            logger.info("Emails have been sent")

    except Exception as exc:
        logger.error(f"Mailing error : {exc}")
        return False

    return True
