
import logging

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from biostar.forum.search import index_posts

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--reindex', action='store_true', default=False, help="Re-index from scratch.")

    def handle(self, *args, **options):

        # Index all posts that have a root.
        posts = Post.objects.exclude(root=None)
        reindex = options['reindex']
        index_posts(posts=posts, reindex=reindex)
