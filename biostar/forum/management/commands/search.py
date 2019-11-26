
import logging
import time
from itertools import count, islice

from django.core.management.base import BaseCommand
from django.conf import settings

from biostar.forum.models import Post
from biostar.forum.search import preform_search, more_like_this

logger = logging.getLogger('engine')


def time_func(func, kwargs):
    """
    Time a function. Return the time it took to finish and subsequent result.
    """

    start = time.time()
    res = func(**kwargs)
    finish = time.time() - start

    return finish, res


class Bunch():
    def __init__(self, **kwargs):
        self.title = ""
        self.content = ""
        self.author_name = ""
        self.author_uid = ""
        self.author_handle = ""
        self.post_type = ""
        self.post_uid = ""
        self.lastedit_date = ""
        self.tags = ""
        self.__dict__.update(kwargs)


def parse_result(result):
    "Return a bunch object for result."

    # Result is a database object.
    if isinstance(result, Post):
        bunched = Bunch(title=result.title, content=result.content,
                        author_name=result.author.profile.name, author_uid=result.author.profile.uid,
                        author_handle=result.author.username, post_type=result.get_type_display(),
                        post_uid=result.uid, lastedit_date=result.lastedit_date,
                        tags=result.tag_val)

    # Result is a dictionary returned from indexed searched
    else:
        bunched = Bunch(title=result['title'], content=result['content'],
                        author_name=result['author'], author_uid=result['author_uid'],
                        author_handle=result['author_handle'], post_type=result['type_display'],
                        post_uid=result['uid'], lastedit_date=result['lastedit_date'],
                        tags=result['tags'])

    return bunched


def printer(result, verbosity=0):
    print(f'Type\t{result.post_type}')
    print(f'Uid\t{result.post_uid}')
    print(f'Title\t{result.title}')
    print(f'Last edited\t{result.lastedit_date}')
    print(f'Author name:\t{result.author_name}')
    print(f'Author handle (username):\t{result.author_name}')
    if verbosity >= 2:
        print(f'Tags\t{result.tags}')
        print(f'Content\t{result.content}')
        print(f'Author name\t{result.author_name}')
        print(f'Author handle (username)\t{result.author_name}')


def print_results(results, limit=10, verbosity=0, finish_time=0.0):

    print('-'*20)
    print(f'Search Time\t{finish_time}')
    print(f'Total results\t{len(results)}')

    stream = zip(count(1), results)
    stream = islice(stream, limit)
    if verbosity <= 0:
        return

    print(f'Showing top {limit} search results.')
    print()
    for index, result in stream:
        result = parse_result(result)
        printer(result=result, verbosity=verbosity)
        print("-" * 2)
        print()

    print('-' * 20)
    return


class Command(BaseCommand):
    help = 'Preform a search and generate report on results.'

    def add_arguments(self, parser):

        parser.add_argument('-db', '--db_search', action='store_true', default=False, help="Preform a database.")
        parser.add_argument('--uid', type=str, required=False,
                            help="Search for posts similar to this one (more list this).")
        parser.add_argument('-q', '--query', type=str, required=False, default='test search index',
                            help="Search for posts similar to this one (more list this).")
        parser.add_argument('-l', '--limit', type=int, default=10,
                            help="Print limited amount of results.")
        parser.add_argument('-vr', '--verbose', type=int, default=0, help="Verbosity level of the report.")

    def handle(self, *args, **options):
        logger.info(f"Database: {settings.DATABASE_NAME}. Index : {settings.INDEX_DIR}")

        db_search = options['db_search']
        uid = options['db_search']
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

            finish, results = time_func(more_like_this, kwargs=dict(uid=uid, db_search=db_search))
            logger.info(f"Post uid: {uid}. Database search: {db_search}")
            print_results(results=results, limit=limit, verbosity=verbosity, finish_time=finish)

            return

        finish, results = time_func(preform_search, kwargs=dict(query=query, db_search=db_search))
        logger.info(f"Query: {query}. Database search: {db_search}")
        print_results(results=results, limit=limit, verbosity=verbosity, finish_time=finish)


