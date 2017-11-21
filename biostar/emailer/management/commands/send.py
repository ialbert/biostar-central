import logging
from django.core.management.base import BaseCommand
from django.contrib.sites.models import Site
from django.contrib.auth import get_user_model
from biostar.emailer import sender
from mailer.engine import send_all

from django.conf import settings


logger = logging.getLogger("engine")



class Command(BaseCommand):
    help = 'Send an email to one person or a whole mailing list (--group_name)'

    def add_arguments(self, parser):
        parser.add_argument('--group_name', type=str, required=False,
                            help="Subscription group name to send a batch email to.")

        parser.add_argument('--from', type=str, required=False,
                            default="mailer@biostars.org", help="The sender email")

        parser.add_argument('--to', type=str, required=False,
                            default="mailer@biostars.org", help="The target email ( for one person)")

        parser.add_argument('--name', type=str, required=False,
                            default="mailer@biostars.org", help="The target name ( for one person)")

        parser.add_argument('--template', type=str, required=False,
                            default="test_email.html", help="Email template sent to mailing-list or target recipient. ")


    def handle(self, *args, **options):
        template_name = options['template']
        group = options['group']
        target_email = options['to']
        sender_email = options["from"]

        # Sender requires a list.
        if not group:

            recipient_list =  target_email.split(",")


        # The object that parsers the template.
        email = sender.EmailTemplate(template_name)

        # This is the context passed to each template.
        site = Site.objects.get_current()
        context = dict(site=site, protocol=settings.PROTOCOL, target_email=target_email)

        # Generate a log message.
        logger.info(f"generating email from:{from_email} to:{target_email} using template:{template_name}")

        # Queues the email into the database.
        email.send(context=context, from_email=sender_email, recipient_list=recipient_list)


        return

