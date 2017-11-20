import logging
from django.core.management.base import BaseCommand
from django.conf import settings


logger = logging.getLogger("engine")



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
        return
