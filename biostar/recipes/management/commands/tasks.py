
from django.core.management.base import BaseCommand
from biostar.utils.helpers import pg_dump
import logging

logger = logging.getLogger('biostar')

DUMP = 'pg_dump'

CHOICES = [ DUMP]


class Command(BaseCommand):
    help = 'Preform action on list of posts.'

    def add_arguments(self, parser):
        parser.add_argument('--action', '-a', type=str, required=True, choices=CHOICES, default='', help='Action to take.')
        parser.add_argument('--user', dest='pg_user', default="www", help='postgres user default=%default')
        parser.add_argument('--prog', dest='prog', default="/usr/local/bin/pg_dump", help='the postgres program default=%default')
        parser.add_argument('--outdir', dest='outdir', default="~/data/", help='output directory default=%default')
        parser.add_argument('--hourly', dest='hourly', action='store_true', default=False, help='hourly datadump'),

    def handle(self, *args, **options):
        action = options['action']

        opts = {DUMP: pg_dump}

        func = opts[action]

        #models.Award.objects.all().delete()

        func(**options)
