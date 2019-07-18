
import logging

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from biostar.forum.search import crawl

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--reindex', action='store_true', default=False, help="Re-index from scratch.")
        parser.add_argument('--overwrite', action='store_true', default=False, help="Overtwrites exisiting index.")
        parser.add_argument('--limit', type=int, default=1000, help="Limits the number of posts")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        reindex = options['reindex']
        overwrite = options['overwrite']
        limit = options['limit']

        # Crawl through posts in batches and index
        crawl(limit=limit, reindex=reindex, overwrite=overwrite)
