from biostar import VERSION
from django.conf import settings
from biostar.forum.models import PostView
from django.core.cache import cache
from datetime import timedelta
from django.db import NotSupportedError

from . import util

def get_traffic(key='traffic', timeout=300, minutes=60):
    """
    Obtains the number of distinct IP numbers.
    """
    traffic = cache.get(key)
    if not traffic:
        recent = util.now() - timedelta(minutes=minutes)
        try:
            traffic = PostView.objects.filter(date__gt=recent).distinct('ip').count()
        except NotSupportedError as exc:
            traffic = PostView.objects.filter(date__gt=recent).values_list('ip')
            traffic = [t[0] for t in traffic]
            traffic = len(set(traffic))
        cache.set(key, traffic, timeout)

    return traffic


def forum(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''

    #print(request.session.get('res'), 'res', request.COOKIES)
    res = request.COOKIES.get('resolution', 'x')
    width, height = res.split('x')

    params = dict(user=request.user, width=width, height=height,
                  TRAFFIC=get_traffic(),
                  VERSION=VERSION, request=request, site_name=settings.SITE_NAME,
                  site_domain=settings.SITE_DOMAIN, google_tracker=settings.GOOGLE_TRACKER)
    return params
