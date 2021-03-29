
from django.conf import settings
from django.db.models import Max, Count
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


def delete_repeats():

    blogs = Blog.objects.annotate(count=Count("blogpost__id"))
    # Order by most blog posts.
    blogs = blogs.order_by("-count")

    seen = set()
    for blg in blogs:
        # Delete the blog if already seen
        if blg.feed in seen:
            blg.delete()
            logger.debug(f"deleted {blg.feed}")
        seen.update([blg.feed])


def dropall():

    Blog.objects.all().delete()

    logger.info("deleted all blogs.")

    return


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--add', dest='add', help='adds blogs to the database')
        parser.add_argument('--download', dest='download', action="store_true", default=False,
                            help='downloads latest feeds')
        parser.add_argument('--report', action='store_true', default=False, help="Reports on the content of the index.")
        parser.add_argument('--update', dest='update', default=0, type=int, help='updates existing blogs with latest feeds')
        parser.add_argument('--fake', dest='fake',  action="store_true", default=False, help='Create fake blog entries.')
        parser.add_argument('--drop', dest='drop',  action="store_true", default=False,
                            help='Delete repeated blogs in database.')
        parser.add_argument('--limit', dest='limit', default=1, type=int, help='object limits')

    def handle(self, *args, **options):
        # Create the planet directory if it is missing
        os.makedirs(settings.PLANET_DIR, exist_ok=True)

        limit = options['limit']
        fname = options['add']
        update = options['update']
        download = options['download']
        drop = options['drop']

        if options['fake']:
            fake()

        fname = options['add']

        if drop:
            dropall()

        if fname:
            auth.add_blogs(fname)

        if download:
            auth.download_blogs()

        if update:
            auth.update_entries(update)
