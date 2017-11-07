import logging
from django.core.management.base import BaseCommand
from django.conf import settings
from django.contrib.sites.models import Site
from django.contrib.auth import get_user_model
from biostar.emailer import sender
from mailer.engine import send_all

#TODO: should it get the engine level logger.
logger = logging.getLogger("biostar")


class Command(BaseCommand):
    help = 'tests email settings'

    def add_arguments(self, parser):
        parser.add_argument('--to', type=str, required=False,
                            default="2@lvh.me", help="The target email")

        parser.add_argument('--from', type=str, required=False,
                            default="mailer@biostars.org", help="The sender email")

        parser.add_argument('--template', type=str, required=False,
                            default="test_email.html", help="The template to use.")

    def handle(self, *args, **options):
        template_name = options['template']
        from_email = options['from']
        target_email = options['to']

        # Sender requires a list.
        recipient_list =  target_email.split(",")

        # The object that parsers the template.
        email = sender.EmailTemplate(template_name)

        # This is the context passed to each template.
        site = Site.objects.get_current()
        context = dict(site=site, protocol=settings.PROTOCOL, target_email=target_email)

        # Generate a log message.
        logger.info(f"generating email from:{from_email} to:{target_email} using template:{template_name}")

        # Queues the email into the database.
        email.send(context=context, from_email=from_email, recipient_list=recipient_list)

        # This sends the accumulated email.
        send_all()
