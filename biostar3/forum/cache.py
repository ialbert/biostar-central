"""
Cache related functionality.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
from django.conf import settings
from django.contrib import messages
from django.core.cache import cache
from .models import *
from datetime import datetime, timedelta

# Group information stored in the cache
GROUP_CACHE = "group-%s"
TRAFFIC_CACHE = "traffic"

from .models import PostView

def get_traffic(request, minutes=60):
    """
    Traffic as post views in the last 60 min.
    """
    global TRAFFIC_CACHE

    count = cache.get(TRAFFIC_CACHE)
    if not count:
        # Set the cache for traffic.
        now = right_now()
        start = now - timedelta(minutes=minutes)
        try:
            count = PostView.objects.filter(date__gt=start).exclude(date__gt=now).distinct(
                'ip').count()
        except NotImplementedError:
            count = PostView.objects.filter(date__gt=start).exclude(date__gt=now).count()

        cache.set(TRAFFIC_CACHE, count, timeout=600)
    return count


def get_group(domain):
    global GROUP_CACHE

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
