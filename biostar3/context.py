from __future__ import absolute_import, division, print_function, unicode_literals

import logging

from django.conf import settings
from biostar3.forum.cache import get_traffic
from biostar3 import VERSION
from biostar3.forum import query
from biostar3.forum.models import Post, PostView, Vote, Message, Award
from django.utils.timezone import utc
from datetime import datetime, timedelta

# The count keys are stored in the cache for the user.
COUNT_KEY_PATT = "count-%s"

# When to refresh user data
USER_SESSION_TIMEOUT = 60 * 10

logger = logging.getLogger('biostar')


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def ago(hours=0, minutes=0, days=0):
    since = now() - timedelta(days=days, hours=hours, minutes=minutes)
    return since


def get_shortcuts(request):
    return settings.DEFAULT_SHORTCUTS


ANON_COUNTS = dict(
    post_count=0, book_count=0,
    vote_count=0, mesg_count=0, badge_count=0,
)

SESSION_COUNT_KEY = "counts"


def get_counts(request, counts=None):
    """
    Returns a dictionary with
    """
    if request.user.is_anonymous():
        return ANON_COUNTS

    sess = request.session

    user = request.user
    count_key = COUNT_KEY_PATT % user.id

    last_login = user.profile.last_login

    # Time has passed.
    elapsed = (now() - last_login).seconds

    # The sessions has not been set yet.
    missing = SESSION_COUNT_KEY not in sess

    # The user session will be updated.
    if elapsed > settings.SESSION_UPDATE_SECONDS or missing:
        counts = dict(
            mesg_count=Message.objects.filter(user=user, unread=True).count(),

            vote_count=Vote.objects.filter(post__author=user, type__in=(Vote.BOOKMARK, Vote.UP)).count(),

            new_vote_count=Vote.objects.filter(post__author=user, type__in=(Vote.BOOKMARK, Vote.UP),
                                               unread=True).count(),

            post_count=Post.objects.filter(author=user).count(),
            book_count=Vote.objects.filter(author=user, type=Vote.BOOKMARK).count(),
            award_count=Award.objects.filter(user=user).count(),
        )
        sess[SESSION_COUNT_KEY] = counts
        user.profile.last_login = now()
        user.profile.save()

    counts = sess.get(SESSION_COUNT_KEY, {})
    return counts


def extras(request):
    # These values will be added to each context

    context = {
        "BIOSTAR_VERSION": VERSION,
        "user": request.user,
        "group": request.group,
        "request": request,
        "recaptcha": settings.RECAPTCHA_PUBLIC_KEY,
        "counts": get_counts(request),
        "shortcuts": get_shortcuts(request),
        "TRAFFIC": get_traffic(request),
        "recent_votes": query.recent_votes(),
        "recent_users": query.recent_users(),
        "recent_awards": query.recent_awards(),
        "recent_replies": query.recent_replies(),
    }

    return context


