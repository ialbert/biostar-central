
import logging
import time
from itertools import count, islice

from django.core.management.base import BaseCommand
from django.conf import settings

from biostar.forum.models import Post
from biostar.forum.search import more_like_this, perform_search

logger = logging.getLogger('engine')


def time_func(func, kwargs):
    """
    Time a function. Return the time it took to finish and subsequent result.
    """

    start = time.time()
    res = func(**kwargs)
    finish = time.time() - start

    return finish, res


def printer(result, verbosity=0):
    print(f'Title\t{result.title}')
    print(f'Type\t{result.type_display}')
    print(f'Uid\t{result.uid}')
    print(f'Last edited\t{result.lastedit_date}')
    print(f'Author name:\t{result.author}')
    print(f'Author handle (username):\t{result.author_handle}')
    if verbosity >= 2:
        print(f'Tags\t{result.tags}')
        print(f'Content\t{result.content}')


def print_results(results, limit=10, verbosity=0, finish_time=0.0, query=''):

    stream = zip(count(1), results)
    stream = islice(stream, limit)
    print('-' * 20)
    if verbosity >= 1:
        print(f'Showing top {limit} search results.')
        print()
        for index, result in stream:
            printer(result=result, verbosity=verbosity)
            print("-" * 2)
            print()
        print(f'Showing top {limit} search results.')

    finish_time = '{:.06f}'.format(finish_time)
    print(f'Search query/uid\t{query}')
    print(f'Search Time\t{finish_time} secs')
    print(f'Total results\t{len(results)}')

    print('-' * 20)
    return


class Command(BaseCommand):
    help = 'Preform a search and generate report on results.'

    def add_arguments(self, parser):

        parser.add_argument('--uid', type=str, required=False,
                            help="Search for posts similar to this one (more list this).")
        parser.add_argument('-q', '--query', type=str, required=False, default='test search index',
                            help="Search for posts similar to this one (more list this).")
        parser.add_argument('-l', '--limit', type=int, default=10,
                            help="Print limited amount of results.")
        parser.add_argument('-vr', '--verbose', type=int, default=1, help="Verbosity level of the report.")

    def handle(self, *args, **options):
        logger.info(f"Database: {settings.DATABASE_NAME}. Index : {settings.INDEX_DIR}")

        uid = options['uid']
        query = options['query']
        limit = options['limit']
        verbosity = options['verbose']

        # Preform a more like this search for a given uid
        if uid:
            post = Post.objects.filter(uid=uid).first()
            if not post:
                logger.info(f"Post uid does not exist : {uid}")
                return

            logger.info(f"Searching for similar posts: {post.title}")
            finish, results = time_func(perform_search, kwargs=dict(query=uid, fields=['uid']))
            logger.info(f"Post uid: {uid}.")
            print_results(results=results, limit=limit, verbosity=verbosity, query=uid, finish_time=finish)

            return

        finish, results = time_func(perform_search, kwargs=dict(query=query))
        logger.info(f"Query: {query}")
        print_results(results=results, limit=limit, query=query, verbosity=verbosity, finish_time=finish)


