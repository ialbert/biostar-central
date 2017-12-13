import logging, os
from django.core.management.base import BaseCommand
from django.contrib.sites.models import Site

from biostar.emailer import sender, models
from mailer.engine import send_all

from django.conf import settings

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)

class Command(BaseCommand):
    help = 'Send an email to a group.'

    def add_arguments(self, parser):
        parser.add_argument('--name', type=str, required=False,
                            help="Subscription group name to send the email to.")

        parser.add_argument('--subject', type=str, required=False, default="email subject",
                            help="Subject of the email.")

        parser.add_argument('--from', type=str, required=False,
                            default="", help="The sender email")

        parser.add_argument('--template', type=str, required=False,
                            default="test_email.html", help="Email template.")


    def handle(self, *args, **options):

        template_name = options['template']
        group_name = options['name']
        subject = options['subject']
        from_email = options['from'] or settings.ADMINS[0][1]

        group = models.EmailGroup.objects.filter(name=group_name).first()
        # Sender requires a list.
        if not group:
            logger.error(f"group.name={group_name} does not exist.")
            return

        # Get the recipients.
        recipients = [person.address.email for person in group.subscription_set.all()]

        logger.info(f"Email list size: {len(recipients)}")

        # Test the templates
        if os.path.isfile(template_name):
            logger.error(f"Missing template: {template_name}")
            return

        # Generate emails.
        logger.info(f"Emails from={from_email} to group.name={group.name} using template:{template_name}")

        # The object that parsers the template.
        email = sender.EmailTemplate(template_name)

        # This is the context passed to each template.
        site = Site.objects.get_current()

        # Accumulate the emails into the database.
        for address in recipients:
            context = dict(site=site, protocol=settings.PROTOCOL, subject=subject)
            email.send(context=context, from_email=from_email, recipient_list=[address])

        # Send the emails.
        send_all()
        logger.info("Emails have been sent")



