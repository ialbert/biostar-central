
import logging
import os
from django.conf import settings
from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from biostar.forum.search import index_posts

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--new', type=bool, default=False, help="Delete existing index and start a fresh one.")

    def handle(self, *args, **options):

        create_new = options['new']

        # Index all posts.
        posts = Post.objects.all()

        #index_dir = os.path.abspath(index_dir)

        # Ensure the index directory exists.
        #os.makedirs(index_dir, exist_ok=True)

        #if not os.path.exists(index_dir):
        #    raise Exception("Index directory does not exist.")

        # Fetch the index

        index_posts(posts=posts)
