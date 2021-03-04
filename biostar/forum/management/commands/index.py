
import logging
from typing import Any

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from django.conf import settings
from biostar.forum import search, spam, models

logger = logging.getLogger('engine')


def clear_spam():
    """
    Clear spam from
    """
    # Get all of the spam posts
    posts = models.Post.objects.filter(spam=models.Post.SPAM).values_list('uid', flat=True)

    nposts = len(posts)

    logger.info(f"{nposts} spam found.")

    for uid in posts:
        spam.remove_spam(uid=uid)

    # Search for each post

    return


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--reset', action='store_true', default=False, help="Resets the indexed flags.")
        parser.add_argument('--remove', action='store_true', default=False, help="Removes the existing index.")
        parser.add_argument('--report', action='store_true', default=False, help="Reports on the content of the index.")
        parser.add_argument('--index', type=int, default=0, help="How many posts to index")
        parser.add_argument('--clear_spam',  action='store_true', default=False, help="Clear search index of spam posts.")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        logger.info(f"Database: {settings.DATABASE_NAME}")
        reset = options['reset']
        remove = options['remove']
        report = options['report']
        index = options['index']
        clear = options['clear_spam']

        # Sets the un-indexed flags to false on all posts.
        if reset:
            logger.info(f"Setting indexed field to false on all post.")
            Post.objects.valid_posts(indexed=True).exclude(root=None).update(indexed=False)

        # Index a limited number yet unindexed posts
        if index:

            # How many total posts can be indexed
            start_count = Post.objects.valid_posts(indexed=False).exclude(root=None).count()
            logger.info(f"Starting with {start_count} unindexed posts")

            posts = Post.objects.valid_posts(indexed=False).exclude(root=None)[:index]
            target_count = len(posts)

            logger.info(f"Indexing {target_count} posts")

            # The list of posts to update
            ids = [ post.id for post in posts ]

            # Add post to search index.
            search.index_posts(posts=posts, overwrite=remove)

            # Set the indexed field to true.
            Post.objects.filter(id__in=ids).update(indexed=True)

            count = Post.objects.valid_posts(indexed=False).exclude(root=None).count()
            logger.info(f"Finished with {count} unindexed posts remaining")

        # Report the contents of the index
        if report:
            search.print_info()

        if clear:
            clear_spam()

