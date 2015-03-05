from __future__ import print_function, unicode_literals, absolute_import, division
import logging

from django.core.management.base import BaseCommand
from django.conf import settings
from django.core.mail import send_mail
from biostar3.forum.mailer import EmailTemplate

logger = logging.getLogger('biostar')

class Command(BaseCommand):
    help = 'Tests email settings'

    def handle(self, *args, **options):

        from_email = settings.DEFAULT_FROM_EMAIL

        # Send email to the first admin user
        name, email = settings.ADMINS[0]

        to = [email]

        data = dict(name="John Doe", content="world")

        em = EmailTemplate("mailer_test.html", data=data)

        recp = ", ".join(to)

        try:
            em.send(from_email=from_email, to=to)
            logger.info("email sent to %s via %s " % (recp, settings.EMAIL_BACKEND))
        except Exception, exc:
            logger.error("email error %s when sending to %s via %s " % (exc, recp, settings.EMAIL_BACKEND))



