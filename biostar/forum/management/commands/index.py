
import logging

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from biostar.forum.search import index_posts

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--reindex', action='store_true', default=False, help="Re-index from scratch.")
        parser.add_argument('--delete', action='store_true', default=False, help="Set posts to un-indexed.")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        posts = Post.objects.exclude(root=None, indexed=False)
        reindex = options['reindex']
        delete = options['delete']

        if delete:
            total = posts.count()
            logger.info(f"Settings {total} posts to un-indexed.")
            posts.update(indexed=False)
            logger.info(f"Un-indexed  {total} posts.")
        else:
            index_posts(posts=posts, reindex=reindex)
