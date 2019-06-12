import logging
from django.core.management.base import BaseCommand

from biostar.emailer import  models, tasks

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

        from_email = options['from'] or settings.ADMINS[0][1]

        group = models.EmailGroup.objects.filter(name=group_name).first()
        # Sender requires a list.
        if not group:
            logger.error(f"group.name={group_name} does not exist.")
            return

        # Get the recipients.
        recipients = [person.address.email for person in group.subscription_set.all()]

        logger.info(f"Email list size: {len(recipients)}")

        tasks.send_email(template_name=template_name, recipient_list=recipients,
                         from_email=from_email)



