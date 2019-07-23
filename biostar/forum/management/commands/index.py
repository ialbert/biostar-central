
import logging

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from biostar.forum import search

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--reset', action='store_true', default=False, help="Resets the indexed flags.")
        parser.add_argument('--remove', action='store_true', default=False, help="Removes the existing index.")
        parser.add_argument('--report', action='store_true', default=False, help="Reports on the content of the index.")
        parser.add_argument('--index', type=int, default=0, help="How many posts to index")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        reset = options['reset']
        remove = options['remove']
        report = options['report']
        index = options['index']

        # Sets the un-indexed flags to false on all posts.
        if reset:
            logger.info(f"Setting indexed field to false on all post.")
            Post.objects.filter(indexed=True).exclude(root=None).update(indexed=False)

        # Index a limited number yet unindexed posts
        if index:
            posts = Post.objects.exclude(root=None, indexed=False)[:index]

            # Add post to search index.
            search.index_posts(posts=posts, overwrite=remove)

            # Set the indexed field to true.
            Post.objects.filter(id__in=posts.values('id')).update(indexed=True)

        # Report the contents of the index
        if report:
            search.print_info()
