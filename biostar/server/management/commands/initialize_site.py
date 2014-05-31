from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
import os, logging
from django.contrib.sites.models import Site
from django.contrib.flatpages.models import FlatPage
from allauth.socialaccount.models import SocialApp, providers

from django.core.exceptions import ImproperlyConfigured
from optparse import make_option

logger = logging.getLogger(__name__)

def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))

class Command(BaseCommand):
    help = 'Initializes content in Biostar'

    def handle(self, *args, **options):
        from biostar import awards
        init_admin()
        init_domain()
        init_social_providers()
        init_flatpages()
        awards.init_awards()

def init_flatpages():
    # list for the flatpages
    names = "faq about help policy api".split()
    site = Site.objects.get_current()
    for name in names:
        url = "/info/%s/" % name
        page = FlatPage.objects.filter(url=url, sites=site)
        if not page:
            path = abspath(settings.FLATPAGE_IMPORT_DIR, name)
            path = "%s.html" % path
            if not os.path.isfile(path):
                logger.error("cannot find flatpage %s" % path)
                continue
            content = file(path).read()
            page = FlatPage.objects.create(url=url, content=content, title=name.capitalize())
            page.sites.add(site)
            page.save()
            logger.info("added flatpage for url: %s" % url)

def init_admin():
    # Add the admin user if it is not present.
    from biostar.apps.users.models import User

    email = settings.ADMIN_EMAIL
    admin = User.objects.filter(id=1)
    if not admin:
        admin = User(
            email=email,
            is_staff=True,
            is_admin=True,
            name=settings.ADMIN_NAME,
            type=User.ADMIN
        )
        admin.set_password(settings.SECRET_KEY)
        admin.save()

        admin.profile.location = settings.ADMIN_LOCATION
        admin.profile.save()

        logger.info(
            "added admin user with email=%s, password=SECRET_KEY, name=%s" % (admin.email, admin.get_full_name()))

def init_domain():
    # Initialize to the current site if it is not present.
    from django.contrib.sites.models import Site

    site = Site.objects.get_current()
    if site.domain != settings.SITE_DOMAIN:
        site.name = settings.SITE_NAME
        site.domain = settings.SITE_DOMAIN
        site.save()
        logger.info("adding site=%s, name=%s, domain=%s" % (site.id, site.name, site.domain))

    # Initialize media folder
    for path in (settings.EXPORT_DIR, settings.MEDIA_ROOT):
        if not os.path.isdir(path):
            os.mkdir(path)


def init_social_providers():
    # Initialize social login providers.

    for name, data in settings.SOCIALACCOUNT_PROVIDERS.items():

        try:
            client_id = data.get('PROVIDER_KEY','')
            secret = data.get('PROVIDER_SECRET_KEY','')
            site = Site.objects.get(id=settings.SITE_ID)

            # Check that the provider is registered
            provider = providers.registry.by_id(name)

            # Code duplication since many2many fields cannot be initialized in one step
            exists = SocialApp.objects.filter(name=name)
            if not exists:
                app = SocialApp(
                    name=name,
                    client_id=client_id,
                    provider=name,
                    secret=secret, key='',
                )
                app.save()
                app.sites.add(site)
                app.save()
                logger.info("initializing social provider %s" % name)

        except Exception, exc:
            raise ImproperlyConfigured("error setting provider %s, %s" % (name, exc))


