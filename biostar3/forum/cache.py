"""
Cache related functionality.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
from django.conf import settings
from django.contrib import messages
from django.core.cache import cache
from .models import *

# Group information stored in the cache
GROUP_CACHE = "group-%s"

def get_group(domain):
    # This is called on every request. Needs to be fast.
    key = GROUP_CACHE % domain
    group = cache.get(key)

    # Found the group in the cache
    if group:
        return group

    # Domain has an alias to default group
    if domain in settings.DEFAULT_SUBDOMAINS:
        domain = settings.DEFAULT_GROUP_DOMAIN

    # Get the group and put it in the cache.
    group = UserGroup.objects.filter(domain__iexact=domain).first()

    # Put the value on into the cache.
    cache.set(key, group, timeout=600)

    return group

# Delete the cache
bust_group_cache = lambda group: cache.delete(GROUP_CACHE % group.domain)
