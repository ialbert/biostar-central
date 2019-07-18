
import logging

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from biostar.forum.search import crawl, init_index

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--reindex', action='store_true', default=False, help="Re-index from scratch.")
        parser.add_argument('--overwrite', action='store_true', default=False, help="Overtwrites exisiting index.")
        parser.add_argument('--limit', type=int, default=1000, help="Limits the number of posts")
        parser.add_argument('-n', '--n_indexed', action='store_true', default=False,
                            help="Return the number of documents already indexed.")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        reindex = options['reindex']
        n_indexed = options['n_indexed']
        overwrite = options['overwrite']
        limit = options['limit']

        if n_indexed:
            ix = init_index()
            length = sum(1 for _ in ix.searcher().documents())
            #length = sum(map(lambda x: 1, ix.searcher().documents()))
            logger.info(f"Documents={length}")
            return

        # Crawl through posts in batches and index
        crawl(limit=limit, reindex=reindex, overwrite=overwrite)
