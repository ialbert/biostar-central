
import logging
from django.db.models import Q
from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from django.conf import settings
from biostar.forum import search, spam, spam2

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--reset', action='store_true', default=False, help="Resets the indexed flags.")
        parser.add_argument('--remove', action='store_true', default=False, help="Removes the existing index.")
        parser.add_argument('--test', action='store_true', default=False,
                            help="Run specificity/sensitivity test against content in database.")
        parser.add_argument('--index', action='store_true', default=False, help="How many posts to index")
        parser.add_argument('--fname', type=str, default='', help="File with spam to index. Separated by --delim.")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        logger.info(f"Database: {settings.DATABASE_NAME}")
        reset = options['reset']
        remove = options['remove']
        index = options['index']
        fname = options['fname']
        test = options['test']

        # Sets the un-indexed flags to false on all posts.
        if reset:
            logger.info(f"Setting indexed field to false on all post.")
            Post.objects.filter(Q(spam=Post.SPAM) | Q(status=Post.DELETED)).update(indexed=False)

        if fname:
            #spam.add_file_to_index(fname=fname, delim=delim)
            logger.info(f"Added {fname} to spam index.")

        # Index a limited number yet unindexed posts
        if index:
            spam.build_spam_index(overwrite=remove)

        # Run specificity and sensitivity tests on posts.
        if test:
            #indexname = "test"
            spam.test_classify()
            spam2.test(indexname=indexname)
            #spam.test_classify()
