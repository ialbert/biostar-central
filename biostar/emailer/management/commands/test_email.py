from __future__ import print_function, unicode_literals, absolute_import, division
import logging

from django.core.management.base import BaseCommand
from django.conf import settings

from biostar.emailer import auth


logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'tests email settings'

    def add_arguments(self, parser):
        parser.add_argument('--emails', type=str,
                            help="Comma separated emails to send test emails to.")
        parser.add_argument('--send', action="store_true",
                            help="Comma separated emails to send test emails to.")

    def handle(self, *args, **options):

        to_emails = options["emails"]
        send = options["send"]
        if to_emails:
            to_emails = to_emails.split(",")

        from_email = settings.DEFAULT_FROM_EMAIL
        subject = "test email "

        recipient_list = to_emails or [settings.ADMINS[0][1]]

        logger.info("sending to %s" % recipient_list)

        auth.notify(template_name="test_email.html", email_list=recipient_list,
                    from_email=from_email, subject=subject, send=send)

