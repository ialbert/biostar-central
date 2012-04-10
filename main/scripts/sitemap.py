"""
Creates a sitemap in the EXPORT directory
"""
import os
from django.conf import settings
from django.contrib.sitemaps import GenericSitemap
from django.contrib.sites.models import Site
from main.server import models, const
from django.utils.encoding import smart_str
from django.template import loader

def path(*args):
    "Generates absolute paths"
    return os.path.abspath(os.path.join(*args))
    
def generate():
    sitemap = GenericSitemap({'queryset': models.Post.objects.filter(type__in=const.POST_TOPLEVEL), })
    urlset = sitemap.get_urls()
    text = loader.render_to_string('sitemap.xml', {'urlset': urlset})
    text = smart_str(text)
    site = Site.objects.get_current()
    fname = path(settings.EXPORT_DIR, 'sitemap.xml')
    print '*** writing sitemap for %s to %s' % (site, fname)
    fp = open(fname, 'wt')
    fp.write(text)
    fp.close()
    print '*** done'
    
if __name__ == '__main__':
    # trigger it on first import        
    generate()