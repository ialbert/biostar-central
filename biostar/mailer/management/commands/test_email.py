import logging
from django.core.management.base import BaseCommand
from django.conf import settings
from django.contrib.sites.models import Site
from django.contrib.auth import get_user_model
from biostar.mailer import sender

logger = logging.getLogger("biostar")

class Command(BaseCommand):
    help = 'tests email settings'

    def add_arguments(self, parser):
        parser.add_argument('--to', type=str, required=False,
                            default="1@lvh.me", help="The target email")
        parser.add_argument('--template', type=str, required=False,
                            default="test_email.html", help="The template to use.")

    def handle(self, *args, **options):

        template_name = options['template']
        target_email= options['to']

        # Get the user model.
        User = get_user_model()

        # Get an admin user. They will be the senders.
        admin = User.objects.filter(is_staff=True).first()

        # The object that parsers the template.
        email = sender.EmailTemplate(template_name)

        # This is the context passed to each template.
        site = Site.objects.get_current()
        context = dict(site=site, user=admin, protocol=settings.PROTOCOL, target_email=target_email)

        # Sending the email.
        logger.info(f"generating email from template {template_name}")

        #
        from_email = admin.email
        # Takes a list of addresss.
        recipient_list = [ target_email ]
        email.send(context=context, from_email=from_email, recipient_list=recipient_list)
