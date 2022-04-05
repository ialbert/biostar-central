import logging, os, csv
from django.core.management.base import BaseCommand
from biostar.emailer import  models, auth
from django.conf import settings


logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)

class Command(BaseCommand):
    help = 'Add an email address to a mailing list'

    def add_arguments(self, parser):
        parser.add_argument('--name', type=str, required=True,
                            help="Subscription group to add emails to.")

        parser.add_argument('--file', type=str, required=True,
                            help="""
                                    Path to text file with email addresses and names. 
                                    Excel csv that contains Name and Email headers.
                                """)

    def handle(self, *args, **options):

        group_name = options['name']
        fname = options["file"]

        group, flag = models.EmailGroup.objects.get_or_create(name=group_name)
        if flag:
            logger.info(f"Created group.name={group}")

        stream = csv.DictReader(open(fname, 'rU'))

        before_count = models.EmailSubscription.objects.all().count()

        for no, row in enumerate(stream):
            email, name = row.get("Email"), row.get("Name")
            if email:
                auth.add_subscription(email=email,name=name, group=group)
            else:
                msg = f'Invalid email at line: {no+1}'
                logger.error(msg)

        after_count = models.Subscription.objects.all().count()

        logger.info(f"Added {after_count-before_count} new subscriptions to group={group_name}")
        logger.info(f"Group={group_name} contains {after_count} subscriptions")
