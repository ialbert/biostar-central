import logging, os
from django.core.management.base import BaseCommand
from biostar.emailer import  models, auth
from django.conf import settings


logger = logging.getLogger("engine")



class Command(BaseCommand):
    help = 'Add an email address to a mailing list'

    def add_arguments(self, parser):
        parser.add_argument('--group', type=str, required=True,
                            help="Subscription group to add emails to.")

        parser.add_argument('--file', type=str, required=True,
                            help="""
                                    Path to text file with email addresses and names. 
                                    Tab-delimited with two columns in each line:
                                    (email_address, name)
                                """)

    def handle(self, *args, **options):

        group = options['group']
        file = options["file"]

        group = models.EmailGroup.objects.filter(name=group).first()
        if not group:
            logger.error(f"Group name {group} does not exist.")
            return

        assert os.path.isfile(file), "String entered is not a file"

        with open(file, "r") as mailing_list:

            for mail in mailing_list:
                email, name = mail.rstrip().split("\t")
                auth.add_sub(email=email, name=name, group=group)

        return

