from __future__ import print_function, unicode_literals, absolute_import, division
# -*- coding: utf8 -*-
__author__ = 'ialbert'
from django.conf import settings
from biostar.apps.accounts import models
from django.db.models.signals import post_syncdb
from django.contrib.sites.models import Site
import logging

logger = logging.getLogger(__name__)

def initialize(sender, **kwargs):

    # Add an admin user if it is not present.
    user, flag = models.User.objects.get_or_create(username="admin")
    if flag:
        user.email = settings.ADMINS[0][1]
        user.is_staff = True
        user.is_superuser = True
        user.first_name = "Señor"
        user.last_name = "Admin"
        user.set_password(settings.SECRET_KEY)
        user.save()
        logger.info("admin username=%s, password=SECRET_KEY, name=%s" % (user.username, user.get_full_name()))

    # Initialize the current site if it is not present.
    site = Site.objects.get_current()
    if site.domain != settings.SITE_DOMAIN:
        site.name = settings.SITE_NAME
        site.domain = settings.SITE_DOMAIN
        site.save()
        logger.info("updated to site=%s, name=%s, domain=%s" % (site.id, site.name, site.domain))

post_syncdb.connect(initialize, sender=models)
