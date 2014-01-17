# -*- coding: utf8 -*-
from __future__ import print_function, unicode_literals, absolute_import, division
from django.conf import settings
from biostar.apps.accounts import models
from django.contrib.sites.models import Site, get_current_site
from allauth.socialaccount.models import SocialApp, providers
from django.core.exceptions import ImproperlyConfigured

import logging

logger = logging.getLogger(__name__)

def initialize(sender, **kwargs):

    # Add the admin user if it is not present.
    email = settings.ADMINS[0][1]
    admin = models.User.objects.filter(email=email)

    if not admin:
        admin = models.User(
            email=email,
            is_staff = True,
            is_admin = True,
            name = "Señor Admin",
        )
        admin.set_password(settings.SECRET_KEY)
        admin.save()
        logger.info("added admin user with email=%s, password=SECRET_KEY, name=%s" % (admin.email, admin.get_full_name()))

    # Initialize the current site if it is not present.
    site = Site.objects.get_current()
    if site.domain != settings.SITE_DOMAIN:
        site.name = settings.SITE_NAME
        site.domain = settings.SITE_DOMAIN
        site.save()
        logger.info("adding site=%s, name=%s, domain=%s" % (site.id, site.name, site.domain))

    # Initialize social login providers.
    for name, data in settings.SOCIALACCOUNT_PROVIDERS.items():

        # not all providers need to have entries
        if "PROVIDER_KEY" not in data:
            continue

        try:
            client_id = data['PROVIDER_KEY']
            secret = data['PROVIDER_SECRET_KEY']
            site = Site.objects.get(id=settings.SITE_ID)
            provider = providers.registry.by_id(name)

            # Code duplication since many2many fields cannot be initialized in one step
            exists = SocialApp.objects.filter(
                name=name, client_id=client_id, provider=name,
                secret=secret, sites=site
            )
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
                logger.info("enabling authentication via %s" % name)

        except Exception, exc:
            raise ImproperlyConfigured("error setting provider %s, %s" % (name, exc))



