"""
Creates a sitemap in the EXPORT directory
"""
import os
from django.conf import settings
from django.contrib.sitemaps import GenericSitemap
from django.contrib.sites.models import Site
from biostar.forum.models import Post
from django.utils.encoding import smart_str
from django.template import loader
from django.core.management.base import BaseCommand, CommandError
import logging
from django.contrib import sitemaps

logger = logging.getLogger("engine")

XML_START = """<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="https://www.sitemaps.org/schemas/sitemap/0.9">"""

XML_END = """
</urlset>"
"""

XML_ROW = """
    <url>
        <loc>https://%s/p/%s/</loc>
    </url>
"""


def ping_google():
    try:
        sitemaps.ping_google('/sitemap.xml')
    except Exception as exc:
        # Bare 'except' because we could get a variety
        # of HTTP-related exceptions.
        logger.error(exc)
        pass


def generate_sitemap():
    site = Site.objects.get_current()

    # Defers all fields beyond uid!
    posts = Post.objects.filter(is_toplevel=True, root__status=Post.OPEN)\
            .exclude(type=Post.BLOG).order_by("-pk").only("uid")

    print(XML_START, end='')
    for post in posts:
        row = XML_ROW % (site.domain, post.uid)
        print(row, end='')
    print(XML_END, end='')


class Command(BaseCommand):
    help = 'Creates a sitemap in the export folder of the site'

    def handle(self, *args, **options):
        generate_sitemap()
        # ping_google()
