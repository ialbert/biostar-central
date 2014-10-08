__author__ = 'ialbert'
from django.conf import settings
from biostar import const, VERSION
from django.core.cache import cache
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post, Vote, PostView
from biostar.apps.badges.models import Award

from datetime import timedelta
from math import pow, e, log
from random import random

LOG2 = log(2)
CACHE_TIMEOUT = settings.CACHE_TIMEOUT


def get_recent_votes():
    votes = Vote.objects.filter(post__status=Post.OPEN).select_related("post").order_by("-date")[
            :settings.RECENT_VOTE_COUNT]
    return votes


def get_recent_users():
    users = User.objects.all().select_related("profile").order_by("-profile__last_login")[:settings.RECENT_USER_COUNT]
    return users


def get_recent_awards():
    awards = Award.objects.all().select_related("user", "badge")
    awards = awards.order_by("-date")[:6]
    return awards


def get_recent_replies():
    posts = Post.objects.filter(type__in=(Post.ANSWER, Post.COMMENT), root__status=Post.OPEN).select_related(("author"))
    posts = posts.order_by("-creation_date")
    posts = posts[:settings.RECENT_POST_COUNT]
    return posts


TRAFFIC_KEY = "traffic"


def get_traffic(minutes=60):
    "Obtains the number of distinct IP numbers "
    global TRAFFIC_KEY
    traffic = cache.get(TRAFFIC_KEY)
    if not traffic:
        recent = const.now() - timedelta(minutes=minutes)
        try:
            traffic = PostView.objects.filter(date__gt=recent).distinct('ip').count()
        except NotImplementedError, exc:
            traffic = PostView.objects.filter(date__gt=recent).values_list('ip')
            traffic = [t[0] for t in traffic]
            traffic = len(set(traffic))
        cache.set(TRAFFIC_KEY, traffic, CACHE_TIMEOUT)
    return traffic


def banner_trigger(request, half=settings.HALF_LIFE):
    user = request.user
    if user.is_anonymous() or user.profile.opt_in:
        return True
    level = pow(e, -LOG2 * user.score/half) + 0.01
    rand = random()
    return rand <= level

def shortcuts(request):
    # These values will be added to each context

    context = {
        "GOOGLE_TRACKER": settings.GOOGLE_TRACKER,
        "GOOGLE_DOMAIN": settings.GOOGLE_DOMAIN,
        "SITE_STYLE_CSS": settings.SITE_STYLE_CSS,
        "SITE_LOGO": settings.SITE_LOGO,
        "SITE_NAME": settings.SITE_NAME,
        "CATEGORIES": settings.CATEGORIES,
        "BIOSTAR_VERSION": VERSION,
        "TRAFFIC": get_traffic(),
        'RECENT_REPLIES': get_recent_replies(),
        'RECENT_VOTES': get_recent_votes(),
        "RECENT_USERS": get_recent_users(),
        "RECENT_AWARDS": get_recent_awards(),
        'USE_COMPRESSOR': settings.USE_COMPRESSOR,
        'COUNTS': request.session.get(settings.SESSION_KEY, {}),
        'SITE_ADMINS': settings.ADMINS,
        'TOP_BANNER': settings.TOP_BANNER,
        'BANNER_TRIGGER': banner_trigger(request),
    }

    return context