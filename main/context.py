"""
Custom context processor
"""
import itertools, random
from datetime import datetime, timedelta

from django.conf import settings

from main import server
from main.server import models
from main.server.const import *
from django.core.cache import cache

TRAFFIC_KEY = 'traffic'
RECENT_TAGS_KEY = "recent-tags"

IMPORTANT_TAGS = models.Tag.objects.filter(name__in=settings.IMPORTANT_TAG_NAMES).order_by('-count')

def get_recent_tags():
    "returns the recent tags"
    posts = models.Post.objects.filter(type=POST_QUESTION, status=POST_OPEN).order_by("-creation_date")[:30]
    tags  = set()
    for p in posts:
        tags.update( p.tag_set.all() )
    tags = list(tags)
    return tags[:30]

def extras(request):
    "Adds more data to each RequestContext"

    user = request.user
        
    # query- this is the primary query that goes into the search box
    q = request.REQUEST.get('q', '')
    
    # match - this is the secondary query value for matching other content
    m = request.REQUEST.get('m', '')
    
    # the tab bar counts
    counts = request.session.get(SESSION_POST_COUNT, {})

    recent_tags = cache.get(RECENT_TAGS_KEY)
    if not recent_tags:
        recent_tags = get_recent_tags()
        cache.set(RECENT_TAGS_KEY, recent_tags, 600)

    # cache the traffic counts
    traffic = cache.get(TRAFFIC_KEY)

    if not traffic:
        try:
            recently = datetime.now() - timedelta(minutes=60)
            traffic = models.PostView.objects.filter(date__gt=recently).distinct('ip').count()
        except Exception, exc:
            traffic = models.PostView.objects.filter(date__gt=recently).count()
        cache.set(TRAFFIC_KEY, traffic, 600)
    
    return { 'BIOSTAR_VERSION': server.VERSION,
             'GOOGLE_TRACKER': settings.GOOGLE_TRACKER,
             'GOOGLE_DOMAIN': settings.GOOGLE_DOMAIN,
             'user':user, 
             'q':q, 
             'm':m,
             'counts':counts,
             'traffic':traffic,
             'important_tags': IMPORTANT_TAGS,
             'recent_tags': recent_tags,
             'params':{}, # this is needed because of the navbar
    }

def popular_tags(request):
    tags1 = models.Tag.objects.filter(name='galaxy') # Special treatment for Galaxy folks
    tags2 = models.Tag.objects.exclude(name='galaxy').order_by('-count')[:5]
    tags = itertools.chain(tags1, tags2)

    return {'popular_tags' : tags}


