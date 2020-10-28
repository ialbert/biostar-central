from __future__ import print_function, unicode_literals, absolute_import, division

import logging

from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.emailer import tasks

logger = logging.getLogger("biostar")

DEFAULT_ADDR = settings.DEFAULT_FROM_EMAIL


class Command(BaseCommand):
    help = 'tests email settings'

    def add_arguments(self, parser):
        parser.add_argument('-f', '--from', type=str, default=DEFAULT_ADDR,
                            help="The sender's email (default=%(default)s")
        parser.add_argument('-t', '--to', type=str, default=DEFAULT_ADDR,
                            help="Comma separated emails of the recipients (default=%(default)s)")

    def handle(self, *args, **options):
        from_email = options["from"]
        recipient_list = options["to"]

        recipient_list = recipient_list.split(",")

        subject = "Test email"

        logger.info(f"settings.EMAIL_BACKEND={settings.EMAIL_BACKEND}")
        logger.info(f"sending test email from {from_email} to {recipient_list}")

        tasks.send_email(template_name="test_email.html", recipient_list=recipient_list,
                         from_email=from_email, subject=subject, name="Testing")

        # Triggers send if the backend is queued.
        tasks.send_all()
