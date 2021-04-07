from datetime import timedelta

from django.conf import settings
from django.core.cache import cache
from django.db import NotSupportedError

from biostar import VERSION
from biostar.accounts.models import is_moderator
from biostar.forum import models
from . import util, const


def get_traffic(key='traffic', timeout=300, minutes=60):
    """
    Obtains the number of distinct IP numbers.
    """
    traffic = cache.get(key)
    if not traffic:
        recent = util.now() - timedelta(minutes=minutes)
        try:
            traffic = models.PostView.objects.filter(date__gt=recent).distinct('ip').count()
        except NotSupportedError as exc:
            traffic = models.PostView.objects.filter(date__gt=recent).values_list('ip')
            traffic = [t[0] for t in traffic]
            traffic = len(set(traffic))
        # It is possible to not have hit any postview yet.
        traffic = traffic or 1
        cache.set(key, traffic, timeout)

    return traffic


def forum(request):
    """
    Additional context applied to each request.
    """

    # Will inject the counts into every session
    if request.user.is_anonymous:
        counts = dict(planet_count=0)
    else:
        counts = request.session.get(settings.SESSION_COUNT_KEY, {})

    params = dict(user=request.user,
                  TRAFFIC=get_traffic(),
                  VERSION=VERSION,
                  request=request,
                  site_name=settings.SITE_NAME,
                  site_domain=settings.SITE_DOMAIN,
                  google_tracker=settings.GOOGLE_TRACKER,
                  IS_MODERATOR=is_moderator(request.user),
                  counts=counts,
                  )

    return params
