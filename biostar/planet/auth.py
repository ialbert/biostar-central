from urllib import request
import os
from datetime import datetime
import logging

import feedparser
from django.utils.encoding import smart_text
from django.utils import timezone
from django.conf import settings

from biostar.planet.models import Blog, BlogPost
from biostar.forum.util import strip_tags, get_uuid

logger = logging.getLogger('engine')


def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))


def create_blogpost(entry, blog):
    date = entry.get('date_parsed') or entry.get('published_parsed')
    date = datetime(date[0], date[1], date[2])
    date = timezone.make_aware(date, timezone=timezone.utc)
    if not entry.title:
        return

    # body = html.clean(entry.description)[:5000]
    body = entry.description
    content = strip_tags(body)
    try:
        post = BlogPost.objects.create(title=entry.title, blog=blog, uid=entry.id, content=content, html=body,
                                       creation_date=date, link=entry.link)

    except Exception as exc:

        logger.error(entry.id)
        logger.warning(f"database error {exc}")
    else:
        logger.info(f"added: {post.title}")

    return


def add_blogpost(blogs, count=3):
    import html
    for blog in blogs:
        logger.debug(f"parsing blog: {blog.id}: {blog.title}")
        try:
            seen = [e.uid for e in BlogPost.objects.filter(blog=blog)]
            seen = set(seen)
            # Parse the blog
            doc = blog.parse()
            # get the new posts
            entries = [e for e in doc.entries if f"{e.id}" not in seen]

            # Only list a few entries
            entries = entries[:count]
            for entry in entries:
                entry.title = smart_text(entry.title)
                entry.title = entry.title.strip()
                # entry.title = html.strip_tags(entry.title)
                entry.title = entry.title.strip()[:200]
                desc = smart_text(entry.description)
                desc = html.unescape(desc)
                entry.description = desc

                # entry.description = html.strip_tags(entry.description)
                create_blogpost(entry=entry, blog=blog)

        except Exception as exc:
            logger.error(f"{exc}")

    return


def blog_from_link(link, user):
    """
    Create a Planet blog post from a Link object.
    """

    return


def add_blog(feed):
    try:
        logger.debug(f"processing {feed}")
        text = request.urlopen(feed).read().decode(errors="surrogatepass")
        doc = feedparser.parse(text)
        title = doc.feed.title
        if hasattr(doc.feed, "description"):
            desc = doc.feed.description
        else:
            desc = ""

        link = doc.feed.link
        blog = Blog.objects.create(title=smart_text(title), feed=feed, link=link, desc=smart_text(desc))

        #add_blogpost(blogs=Blog.objects.filter(id=blog.id))
        blog.download()
        logger.info(f"adding {blog.title}")
        logger.debug(f"link: {blog.link}")
        logger.debug(blog.desc)

    except Exception as exc:
        logger.error(f"error {exc} parsing {feed}")
        blog = None

    return blog


def download_blogs():
    blogs = Blog.objects.filter(active=True, remote=True)
    for blog in blogs:
        logger.info(f"downloading: {blog.title}")
        blog.download()


def update_entries(count=3):
    blogs = Blog.objects.filter(active=True, remote=True)

    # Update blog posts for active blogs
    add_blogpost(blogs=blogs, count=count)


def add_blogs(add_fname):
    # Strip newlines
    urls = map(str.strip, open(add_fname))

    # Attempt to populate the database with the feeds
    for feed in urls:

        exists = Blog.objects.filter(feed=feed).first()

        if exists:
            logger.info(f"Blog already exists feed={feed}")
        else:
            add_blog(feed)
