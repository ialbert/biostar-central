
import logging

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from biostar.forum.search import index_posts

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--reindex', action='store_true', default=False, help="Re-index from scratch.")
        parser.add_argument('--overwrite', action='store_true', default=False, help="Overtwrites exisiting index.")
        parser.add_argument('--limit', type=int, default=1000, help="Limits the number of posts")

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        reindex = options['reindex']
        overwrite = options['overwrite']
        limit = options['limit']

        if reindex:
            logger.info(f"Setting indexed field to false on all post.")
            Post.objects.exclude(root=None, indexed=False).update(indexed=False)

        if overwrite:
            logger.info("Overwriting the index")

        # Index a limited number of posts
        posts = Post.objects.exclude(root=None, indexed=False)[:limit]

        # Add post to search index.
        index_posts(posts=posts, overwrite=overwrite)

        # Set the indexed field to true.
        Post.objects.filter(id__in=posts.values('id')).update(indexed=True)
