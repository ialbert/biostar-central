import logging
from django.core.management.base import BaseCommand
from biostar.emailer import  models
from django.conf import settings


logger = logging.getLogger("engine")




def subscribe(address, group, template=None):


    pass




class Command(BaseCommand):
    help = 'Add an email address to a mailing list'

    def add_arguments(self, parser):
        parser.add_argument('--group', type=str, required=True,
                            help="Subscription list to add emails to.")

        parser.add_argument('--file', type=str, required=True,
                            help="""
                                    Path to text file with email addresses and names. 
                                    Tab-delimited with two columns in each line:
                                    (email_address, name)
                                """)

    def handle(self, *args, **options):
        group = options['group']
        file = options["file"]

        group = ""
        return

