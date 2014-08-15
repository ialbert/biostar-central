__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError
import os, logging
from optparse import make_option

logger = logging.getLogger(__name__)

os.environ['DJANGO_COLORS'] = 'nocolor'

class Command(BaseCommand):
    help = 'Performs actions on users'

    option_list = BaseCommand.option_list + (
        make_option('--award', dest='award', action='store_true', default=False,
                    help='goes over the users and attempts to create awards'),
    )

    def handle(self, *args, **options):

        if options['award']:
            crawl_awards()

def crawl_awards():
    from biostar.apps.users.models import User
    from biostar.awards import create_user_award
    import random

    ids = [ u[0] for u in User.objects.all().values_list("id") ]

    random.shuffle(ids)

    ids = ids[:100]

    for pk  in ids:
        user = User.objects.get(pk=pk)
        #logger.info("%s: %s" % (user.id, user.name))
        create_user_award.delay(user=user)


