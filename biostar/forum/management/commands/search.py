
import logging

from django.core.management.base import BaseCommand
from django.conf import settings
from biostar.forum.models import Post
from biostar.forum.search import preform_search, more_like_this

logger = logging.getLogger('engine')


def print_results(results, limit=10):

    return


class Command(BaseCommand):
    help = 'Preform a search and generate report on results.'

    def add_arguments(self, parser):

        parser.add_argument('--db_search', action='store_true', default=False, help="Preform a database.")
        parser.add_argument('--uid', type=str, required=False,
                            help="Search for posts similar to this one (more list this).")
        parser.add_argument('--query', type=str, required=False, default='test search index',
                            help="Search for posts similar to this one (more list this).")
        parser.add_argument('--limit', type=int, default=10, help="Print limited amount of results.")

    def handle(self, *args, **options):

        logger.info(f"Searching for similar posts to : {settings}")
        db_search = options['db_search']
        uid = options['db_search']
        query = options['query']

        if uid:
            # Preform a more like this search for a given uid
            post = Post.objects.filter(uid=uid).first()
            if not post:
                logger.info(f"Post uid does not exist : {uid}")
                return

            logger.info(f"Searching for similar posts to : {post.title}")
            results = more_like_this(uid=uid, db_search=db_search)

            pass

        results = preform_search(query=query, db_search=db_search)

        print(db_search, uid, query)


