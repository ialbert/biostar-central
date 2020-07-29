import logging
from datetime import datetime
from django.core.management.base import BaseCommand
from django.conf import settings
from biostar.forum.models import Post
from biostar.forum import sync

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Sync contents of one site to this one.'

    def add_arguments(self, parser):

        parser.add_argument('--host', type=str, default="", help="Remote host name housing the postgres database.")
        parser.add_argument('--dbname', type=str, default="", help="Remote postgres database name.")
        parser.add_argument('--port', type=str, default="5432", help="Postgres port on remote host.")
        parser.add_argument('--user', type=str, default="www", help="Postgres user.")
        parser.add_argument('--password', type=str, default="", help="Postgres password.")
        parser.add_argument('--batch', type=int, default=10, help="How many posts to load for the given date range.")
        parser.add_argument('--start', type=str, default="",
                            help="""Start syncing from this date; ISO format. <year>-<month>-<date> eg: 2013-02-14""")
        parser.add_argument('--update', action='store_true', default=True,
                            help="""Preform update when syncing.""")
        parser.add_argument('--reset', action='store_true', default=False,
                            help="""Reset the last_synced date stored in the local database.
                                    This only matters when syncing older dates ( --range is negative ). """)
        parser.add_argument('--report', action='store_true', default=False,
                            help="""Print report on what has/has not been synced between databases. """)

        parser.add_argument('--range', type=int, default=1,
                            help="""Number of days to sync relative to --start.  
                                    If --start is empty and --range is negative ( syncing older posts ) 
                                    last_synced from db is chosen as the start.
                                    
                                    Use negative numbers to indicate looking syncing older posts. 
                                    eg. -1 --> load post 1 day before start date
                                         1 --> load posts 1 day after start date. """)

    def handle(self, *args, **options):
        # Print database that will be synced.
        logger.info(f"Local Database: {settings.DATABASE_NAME}")

        start_date = options['start']
        date_range = options['range']

        report = options['report']

        start_date = datetime.fromisoformat(start_date) if start_date else ''

        if report:
            sync.report(start=start_date, days=date_range, options=options)
        else:
            sync.sync_db(start=start_date, days=date_range, options=options)