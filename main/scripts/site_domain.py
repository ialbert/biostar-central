"""
Update the site settings
"""
from django.conf import settings
from django.contrib.sites.models import Site

def update_site_domain():
    site = Site.objects.get(id=settings.SITE_ID)
    print "*** current site domain %s" % site.domain
    if site.domain != settings.SITE_DOMAIN:
        print '--- setting site domain to %s' % settings.SITE_DOMAIN
        site.domain = settings.SITE_DOMAIN
        site.save()

if __name__ == '__main__':
    # trigger it on first import        
    update_site_domain()