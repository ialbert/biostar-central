
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


def dropall():

    Blog.objects.all().delete()

    logger.info("deleted all blogs.")

    return

def init_local(update):
    """
    Creates a local blog
    """
    blog, created = Blog.objects.get_or_create(title="Local blog", remote=False)
    for step in range(update):
        logger.info("adding local blog post")
        BlogPost.objects.create(blog=blog, title='Local blog post', content="Lorem ipsum", creation_date=now())

class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):

        parser.add_argument('--add', dest='add', help='adds blogs to the database')
        parser.add_argument('--download', dest='download', action="store_true", default=False,
                            help='downloads latest feeds')
        parser.add_argument('--report', action='store_true', default=False, help="Reports on the content of the index.")
        parser.add_argument('--update', dest='update', default=0, type=int, help='updates existing blogs with latest feeds')
        parser.add_argument('--local', dest='local',  action="store_true", default=False, help='Creates a local blog')
        parser.add_argument('--drop', dest='drop',  action="store_true", default=False,
                            help='Delete repeated blogs in database.')


    def handle(self, *args, **options):
        # Create the planet directory if it is missing
        os.makedirs(settings.PLANET_DIR, exist_ok=True)

        fname = options['add']
        update = options['update']
        download = options['download']
        drop = options['drop']
        local = options['local']

        if local:
            init_local(update=update)
            return

        if drop:
            dropall()

        if fname:
            auth.add_blogs(fname)

        if download:
            auth.download_blogs()

        if update:
            auth.update_entries(update)
