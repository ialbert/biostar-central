from __future__ import print_function, unicode_literals, absolute_import, division
import logging

from django.core.management.base import BaseCommand
from django.conf import settings
from django.core.mail import send_mail


logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'tests email settings'

    def handle(self, *args, **options):
        from_email = settings.DEFAULT_FROM_EMAIL
        subject = "[biostar] test email "

        recipient_list = [settings.ADMIN_EMAIL]

        params = dict(subject=subject, from_email=from_email, recipient_list=recipient_list)

        message = """
        Hello,

        this is an email sent via the

        test_email

        Biostar management command. Parameters:

        from_email = %(from_email)s
        recipient_list = %(recipient_list)s
        subject = %(subject)s

        """ % params

        logger.info("sending to %s" % recipient_list)
        send_mail(subject=subject, message=message, from_email=from_email, recipient_list=recipient_list)

