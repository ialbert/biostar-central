
import logging
from typing import Any
import os, sys
from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from django.conf import settings
from biostar.forum import search
from biostar.utils.decorators import check_lock

logger = logging.getLogger('engine')

LOCK = os.path.join(settings.INDEX_DIR, 'flag')


@check_lock(LOCK)
def build(size, remove=False):
    """
    Builds search index
    """

    # Get top level posts that have not been indexed.
    posts = Post.objects.valid_posts(indexed=False, is_toplevel=True).exclude(root=None)[:size]
    target_count = len(posts)

    # The list of posts to update
    ids = [post.id for post in posts]

    # Add post to search index.
    search.index_posts(posts=posts, overwrite=remove)

    # Set the indexed field to true.
    Post.objects.filter(id__in=ids).update(indexed=True)

    count = Post.objects.valid_posts(indexed=False, is_toplevel=True).exclude(root=None).count()

    logger.info(f"Indexed {target_count} posts, {count} unindexed posts remaining")

    # Take spam posts that have been indexed and remove.
    spam_posts = Post.objects.filter(spam=Post.SPAM, indexed=False)[:size]
    sids = [post.id for post in spam_posts]

    for post in spam_posts:
        # Remove spam from search index.
        search.remove_post(post=post)

    # Update the spam indexed flag.
    Post.objects.filter(id__in=sids).update(indexed=True)

    # Add to spam index
    logger.info(f"Removed {len(sids)} spam posts from index")


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--reset', action='store_true', default=False, help="Resets the indexed flags.")
        parser.add_argument('--remove', action='store_true', default=False, help="Removes the existing index.")
        parser.add_argument('--report', action='store_true', default=False, help="Reports on the content of the index.")
        parser.add_argument('--size', type=int, default=0, help="How many posts to index")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        logger.debug(f"Database: {settings.DATABASE_NAME}")
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
            build(size=size, remove=remove)

        # Report the contents of the index
        if report:
            search.print_info()

