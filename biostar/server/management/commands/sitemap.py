"""
Creates a sitemap in the EXPORT directory
"""
import os
from django.conf import settings
from django.contrib.sitemaps import GenericSitemap
from django.contrib.sites.models import Site
from biostar.apps.posts.models import Post
from django.utils.encoding import smart_str
from django.template import loader
from django.core.management.base import BaseCommand, CommandError
import logging
from django.contrib import sitemaps

logger = logging.getLogger("command")

class Command(BaseCommand):
    help = 'Creates a sitemap in the export folder of the site'

    def handle(self, *args, **options):
        generate_sitemap()
        #ping_google()

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))

def ping_google():
    try:
        sitemaps.ping_google('/sitemap.xml')
    except Exception, exc:
        # Bare 'except' because we could get a variety
        # of HTTP-related exceptions.
        logger.error(exc)
        pass

def generate_sitemap():
    sitemap = GenericSitemap({
        'queryset': Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(type=Post.BLOG),
    })
    urlset = sitemap.get_urls()
    text = loader.render_to_string('sitemap.xml', {'urlset': urlset})
    text = smart_str(text)
    site = Site.objects.get_current()
    fname = path(settings.STATIC_ROOT, 'sitemap.xml')
    logger.info('*** writing sitemap for %s to %s' % (site, fname))
    fp = open(fname, 'wt')
    fp.write(text)
    fp.close()
    logger.info('*** done')

if __name__ == '__main__':
    #generate_sitemap()
    pass