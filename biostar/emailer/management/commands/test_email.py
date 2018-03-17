from __future__ import print_function, unicode_literals, absolute_import, division
import logging

from django.core.management.base import BaseCommand
from django.conf import settings

from biostar.emailer import auth


logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'tests email settings'

    def handle(self, *args, **options):
        from_email = settings.DEFAULT_FROM_EMAIL
        subject = "[biostar] test email "

        recipient_list = [ settings.ADMINS[0][1] ]

        logger.info("sending to %s" % recipient_list)

        auth.notify(template_name="test_email.html", email_list=recipient_list,
                    from_email=from_email, subject=subject, send=False)

