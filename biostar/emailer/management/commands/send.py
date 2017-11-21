import logging
from django.core.management.base import BaseCommand
from django.contrib.sites.models import Site

from biostar.emailer import sender, models
from mailer.engine import send_all

from django.conf import settings


logger = logging.getLogger("engine")



class Command(BaseCommand):
    help = 'Send an email to one person or a whole mailing list (--group_name)'

    def add_arguments(self, parser):
        parser.add_argument('--group', type=str, required=False,
                            help="Subscription group name to send a batch email to.")

        parser.add_argument('--from', type=str, required=False,
                            default="mailer@biostars.org", help="The sender email")

        parser.add_argument('--to', type=str, required=False,
                            default="mailer@biostars.org", help="The target email")

        parser.add_argument('--name', type=str, required=False,
                            default="mailer", help="The target name")

        parser.add_argument('--subject', type=str, required=False,
                            default="mailer", help="Email subject")

        parser.add_argument('--template', type=str, required=False,
                            default="base_email.html", help="Email template sent to mailing-list or target recipient. ")


    def handle(self, *args, **options):

        template_name = options['template']
        group = options['group']
        target_email = options['to']
        sender_email = options["from"]
        subject = options["subject"]

        group = models.EmailGroup.objects.filter(name=group).first()

        # Sender requires a list.
        if not group:
            recipient_list =  target_email.split(",")
        else:
            recipient_list = [person.address.email for person in group.subscription_set.all()]

        # The object that parsers the template.
        email = sender.EmailTemplate(template_name)

        # This is the context passed to each template.
        site = Site.objects.get_current()
        context = dict(site=site, protocol=settings.PROTOCOL, target_email=target_email,
                       subject=subject, group=group)

        logger.info(f"generating emails from:{sender_email} to:{group or target_email} using template:{template_name}")

        # Queues the emails into the database.
        email.send(context=context, from_email=sender_email, recipient_list=recipient_list)

        #TODO: send_all doesnt really work; need to send_mass_mail maybe?
        # This sends the accumulated email.
        send_all()


