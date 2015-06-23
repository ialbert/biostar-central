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


