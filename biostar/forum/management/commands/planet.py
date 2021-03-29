
from django.conf import settings

from biostar.forum.util import now
from biostar.planet.models import Blog, BlogPost
from django.core.management.base import BaseCommand
import os
import logging

from biostar.planet import auth

logger = logging.getLogger('engine')


def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))


def fake():
    """
    Create a fake blog post.
    """

    # Get or create a blog
    blog, created = Blog.objects.get_or_create(title="Fake")

    BlogPost.objects.create(blog=blog, title='Creating a fake blog post.', creation_date=now())

class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--add', dest='add', help='adds blogs to the database')
        parser.add_argument('--download', dest='download', action="store_true", default=False,
                            help='downloads latest feeds')
        parser.add_argument('--report', action='store_true', default=False, help="Reports on the content of the index.")
        parser.add_argument('--update', dest='update', default=0, type=int, help='updates existing blogs with latest feeds')
        parser.add_argument('--fake', dest='fake',  action="store_true", default=False, help='Create fake blog entries.')

    def handle(self, *args, **options):
        # Create the planet directory if it is missing
        os.makedirs(settings.PLANET_DIR, exist_ok=True)

        if options['fake']:
            fake()

        fname = options['add']
        if fname:
            fname = os.path.abspath(os.path.join(settings.PLANET_DIR, fname))
            auth.add_blogs(fname)

        if options['download']:
            fname = os.path.abspath(os.path.join(settings.PLANET_DIR, 'example-feeds.txt'))
            auth.add_blogs(fname)
            auth.download_blogs()

        count = options['update']
        if count:
            auth.update_entries(count)
