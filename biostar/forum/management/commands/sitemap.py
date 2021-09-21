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

URLSET_START = """<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="https://www.sitemaps.org/schemas/sitemap/0.9">
"""

URLSET_END = """
</urlset>
"""

URLSET_ROW = """
    <url>
        <loc>https://%s/p/%s/</loc>
        <lastmod>%s</lastmod>
    </url>
"""

SITEMAP_XML = """<?xml version="1.0" encoding="UTF-8"?>
<sitemapindex xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
%s
</sitemapindex>
"""

SITEMAP_ROW = """
    <sitemap>
        <loc>https://%s/static/sitemap_%d.xml</loc>
    </sitemap>
"""


def ping_google():
    try:
        sitemaps.ping_google('/sitemap.xml')
    except Exception as exc:
        # Bare 'except' because we could get a variety
        # of HTTP-related exceptions.
        logger.error(exc)
        pass


def generate_sitemap(index, batch):
    site = Site.objects.get_current()

    # Generates the sitemap index.
    if index:
        body = []
        for step in range(index):
            body.append(SITEMAP_ROW % (site.domain, step+1))
        text = "".join(body)
        print(SITEMAP_XML % text)

        return

    if batch:
        N = 50000
        # Defers all fields beyond uid!
        start = (batch-1) * N
        end = batch * N
        posts = Post.objects.valid_posts(is_toplevel=True, root__status=Post.OPEN) \
            .exclude(type=Post.BLOG).order_by("-pk").only("uid", "lastedit_date")[start:end]

        print(URLSET_START, end='')
        for post in posts:
            lastmod = post.lastedit_date.strftime("%Y-%m-%d")
            row = URLSET_ROW % (site.domain, post.uid, lastmod)
            print(row, end='')
        print(URLSET_END, end='')


class Command(BaseCommand):
    help = 'Creates a sitemap in the export folder of the site'

    def add_arguments(self, parser):
        parser.add_argument('--index', default=0, help="Writes an index")
        parser.add_argument('--batch', default=0, help="50K URL in a batch")

    def handle(self, *args, **options):
        index = int(options['index'])
        batch = int(options['batch'])
        generate_sitemap(index=index, batch=batch)
        # ping_google()
