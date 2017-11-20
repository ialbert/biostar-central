import logging
from django.core.management.base import BaseCommand
from django.conf import settings


logger = logging.getLogger("engine")




def send(address, group):

    return




class Command(BaseCommand):
    help = 'Send an email to one person or a whole mailing list (--group_name)'

    def add_arguments(self, parser):
        parser.add_argument('--group_name', type=str, required=False,
                            default="2@lvh.me", help="Subscription group to send a batch email too")

        parser.add_argument('--from', type=str, required=False,
                            default="mailer@biostars.org", help="The sender email")

        parser.add_argument('--to', type=str, required=False,
                            default="mailer@biostars.org", help="The target email")

        parser.add_argument('--name', type=str, required=False,
                            default="mailer@biostars.org", help="The target name")

        parser.add_argument('--template', type=str, required=False,
                            default="test_email.html", help="Email template that overrides one  ")


    def handle(self, *args, **options):
        template_name = options['template']
        group = options['from']
        target_email = options['to']
        return

