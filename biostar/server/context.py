__author__ = 'ialbert'
from django.conf import settings
from biostar import VERSION
from django.core.cache import cache
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post, Vote, Subscription

CACHE_TIMEOUT = settings.CACHE_TIMEOUT
RECENT_VOTES_KEY = "RECENT_VOTES_KEY"
RECENT_USERS_KEY = "RECENT_USERS_KEY"
RECENT_REPLIES_KEY = "RECENT_REPLIES_KEY"

def get_recent_votes():
    votes = cache.get(RECENT_VOTES_KEY)
    if not votes:
        votes = Vote.objects.filter(post__status=Post.OPEN).select_related("post").order_by("-date")[:settings.RECENT_VOTE_COUNT]
        cache.set(RECENT_VOTES_KEY, votes, CACHE_TIMEOUT)
    return votes

def get_recent_users():
    users = cache.get(RECENT_USERS_KEY)
    if not users:
        users = User.objects.all().order_by("-profile__last_login")[:settings.RECENT_USER_COUNT]
        cache.set(RECENT_USERS_KEY, users, CACHE_TIMEOUT)
    return users

def get_recent_replies():
    posts = cache.get(RECENT_REPLIES_KEY)
    if not posts:
        posts = Post.objects.filter(type__in=(Post.ANSWER, Post.COMMENT))\
                    .select_related("author").order_by("-creation_date")[:settings.RECENT_POST_COUNT]
        cache.set(RECENT_REPLIES_KEY, posts, CACHE_TIMEOUT)
    return posts

def get_traffic():
    return 1

def shortcuts(request):
    # These values will be added to each context

    context = {
        "GOOGLE_TRACKER": settings.GOOGLE_TRACKER,
        "SITE_STYLE_CSS": settings.SITE_STYLE_CSS,
        "SITE_LOGO": settings.SITE_LOGO,
        "CATEGORIES": settings.CATEGORIES,
        "BIOSTAR_VERSION": VERSION,
        "TRAFFIC": get_traffic(),
        'RECENT_REPLIES': get_recent_replies(),
        'RECENT_VOTES': get_recent_votes(),
        'USE_COMPRESSOR': settings.USE_COMPRESSOR,
        'COUNTS': request.session[settings.SESSION_KEY],
        'SITE_NAME': settings.SITE_NAME,
    }

    return context