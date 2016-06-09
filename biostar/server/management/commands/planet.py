__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError
import os, logging
from optparse import make_option
from string import strip
import feedparser, urllib
from django.utils.encoding import smart_text
from datetime import datetime
from django.utils import timezone

logger = logging.getLogger('simple-logger')

def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))


class Command(BaseCommand):
    help = 'Performs planet based data collection'

    option_list = BaseCommand.option_list + (
        make_option('--add', dest='add',
                    help='adds blogs to the database'),
        make_option('--download', dest='download', action="store_true", default=False,
                    help='downloads latest feeds'),
        make_option('--update', dest='update', default=0, type=int,
                    help='updates with latest feeds'),
    )

    def handle(self, *args, **options):
        # Create the planet directory if it is missing
        if not os.path.isdir(settings.PLANET_DIR):
            logger.info("creating planet directory %s" % settings.PLANET_DIR)
            os.mkdir(settings.PLANET_DIR)

        add_fname = options['add']
        if add_fname:
            add_blogs(add_fname)

        if options['download']:
            download_blogs()

        count = options['update']
        if count:
            update_entries(count)

def add_blog(feed):
    from biostar.apps.planet.models import Blog
    # makes it easier to troubleshoot when thing fail
    fname = abspath(settings.PLANET_DIR, 'add-blog.xml')
    try:
        text = urllib.urlopen(feed).read()
        stream = file(fname, 'wt')
        stream.write(text)
        stream.close()
        doc = feedparser.parse(fname)
        title = doc.feed.title
        if hasattr(doc.feed, "description"):
            desc = doc.feed.description
        else:
            desc = ""
        link = doc.feed.link
        blog = Blog.objects.create(title=smart_text(title), feed=feed, link=link, desc=smart_text(desc))
        logger.info("adding %s" % blog.title)
        logger.info("link: %s" % blog.link)
        logger.info(blog.desc)
    except Exception as exc:
        logger.error("error %s parsing %s" % (exc, feed))
        blog = None

    logger.info('-' * 10)
    return blog

def download_blogs():
    from biostar.apps.planet.models import Blog

    blogs = Blog.objects.filter(active=True)
    for blog in blogs:
        logger.info("downloading: %s" % blog.title)
        blog.download()

def update_entries(count=3):
    from biostar.apps.planet.models import Blog, BlogPost
    from biostar.apps.util import html

    #BlogPost.objects.all().delete()

    blogs = Blog.objects.filter(active=True)

    for blog in blogs:
        logger.info("parsing: %s: %s" % (blog.id, blog.title))
        try:
            seen = [e.uid for e in BlogPost.objects.filter(blog=blog)]
            seen = set(seen)

            # Parse the blog
            doc = blog.parse()

            # get the new posts
            entries = [ e for e in doc.entries if e.id not in seen ]

            # Only list a few entries
            entries = entries[:count]

            for r in entries:
                r.title = smart_text(r.title)
                r.title = r.title.strip()
                r.title = html.strip_tags(r.title)
                r.title = r.title.strip()[:200]
                r.description = smart_text(r.description)
                r.description = html.strip_tags(r.description)

                date = r.get('date_parsed') or r.get('published_parsed')
                date = datetime(date[0], date[1], date[2])
                date = timezone.make_aware(date, timezone=timezone.utc)
                if not r.title:
                    continue
                body = html.clean(r.description)[:5000]
                content = html.strip_tags(body)
                try:
                    post = BlogPost.objects.create(title=r.title, blog=blog, uid=r.id, content=content, html=body, creation_date=date, link=r.link)
                except Exception as exc:
                    logger.error(r.title)
                    logger.error("database error %s" % exc)
                else:
                    logger.info("added: %s" % post.title)

        except KeyError as exc:
            logger.error("%s" % exc)


def add_blogs(add_fname):
    from biostar.apps.planet.models import Blog

    #Blog.objects.all().delete()

    # Strip newlines
    urls = map(strip, open(add_fname))

    # Keep feeds with urls that do not yet exists
    urls = filter(lambda url: not Blog.objects.filter(feed=url), urls)

    #urls = urls[:1]

    # Attempt to populate the database with the feeds
    for feed in urls:
        logger.info("parsing %s" % feed)
        add_blog(feed)

