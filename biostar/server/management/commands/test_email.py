from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
import os
from django.core.mail import send_mail
from biostar.apps.util import mailer

class Command(BaseCommand):
    help = 'tests email settings'

    def handle(self, *args, **options):
        from_email = settings.DEFAULT_FROM_EMAIL
        subject = "[biostar] test email"
        message = """
        Test Email.

        Sent via the test_email management command.
        """
        recipient_list = ["istvan.albert@gmail.com"]

        send_mail(subject=subject, message=message, from_email=from_email, recipient_list=recipient_list)

