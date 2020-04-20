
import logging
from django.db.models import Q
from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from django.conf import settings
from biostar.forum import search, spam, api

logger = logging.getLogger('engine')


def load_from_api(base_url, api_key):

    def load_posts():
        return

    def load_users():
        return

    return


class Command(BaseCommand):
    help = 'Sync contents of one site to this one.'

    def add_arguments(self, parser):

        parser.add_argument('--remote_url', type=str, default="", help="Remote url used for api access.")
        parser.add_argument('--api_key', type=str, default="", help="Remote url api key to decode emails.")
        parser.add_argument('--batch', type=str, default="", help="Remote url api key to decode emails.")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        logger.info(f"Database: {settings.DATABASE_NAME}")

        remote_url = options['remote_url']





