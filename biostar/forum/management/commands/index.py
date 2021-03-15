
import logging
from typing import Any

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from django.conf import settings
from biostar.forum import search, spam

logger = logging.getLogger('engine')


def handle_spam(posts):
    """
    Remove list of spam from search index.
    Add posts to spam index.
    """

    for post in posts:
        # Remove spam from search index.
        search.remove_post(post=post)
        # Add spam to its own index
        spam.add_spam(post=post)


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--reset', action='store_true', default=False, help="Resets the indexed flags.")
        parser.add_argument('--remove', action='store_true', default=False, help="Removes the existing index.")
        parser.add_argument('--report', action='store_true', default=False, help="Reports on the content of the index.")
        parser.add_argument('--size', type=int, default=0, help="How many posts to index")
        #parser.add_argument('--clear_spam',  action='store_true', default=False, help="Clear search index of spam posts.")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        logger.info(f"Database: {settings.DATABASE_NAME}")
        reset = options['reset']
        remove = options['remove']
        report = options['report']
        size = options['size']

        # Sets the un-indexed flags to false on all posts.
        if reset:
            logger.info(f"Setting indexed field to false on all post.")
            Post.objects.valid_posts(indexed=True).exclude(root=None).update(indexed=False)

        # Index a limited number yet unindexed posts
        if size:

            posts = Post.objects.valid_posts(indexed=False).exclude(root=None)[:size]
            target_count = len(posts)

            # The list of posts to update
            ids = [ post.id for post in posts ]

            # Add post to search index.
            search.index_posts(posts=posts, overwrite=remove)

            # Set the indexed field to true.
            Post.objects.filter(id__in=ids).update(indexed=True)

            count = Post.objects.valid_posts(indexed=False).exclude(root=None).count()

            logger.info(f"Indexed {target_count} posts, {count} unindexed posts remaining")

            # Take spam posts that have been indexed and remove.
            spam_posts = Post.objects.filter(spam=Post.SPAM, indexed=True)[:size]

            # Remove spam post from index, if present.
            handle_spam(posts=spam_posts)

            # Add to spam index
            logger.info(f"Removed {spam.count()} spam posts, from index")

        # Report the contents of the index
        if report:
            search.print_info()

