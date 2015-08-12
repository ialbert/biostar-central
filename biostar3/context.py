from __future__ import absolute_import, division, print_function, unicode_literals

import logging

from django.conf import settings
from biostar3.forum.cache import get_traffic
from biostar3 import VERSION
from biostar3.forum import query
from biostar3.utils import utils
from biostar3.forum.models import Post, PostView, Vote, Message, Award
from django.utils.timezone import utc
from datetime import datetime, timedelta

# The count keys are stored in the cache for the user.
COUNT_KEY_PATT = "count-%s"

# When to refresh user data
USER_SESSION_TIMEOUT = 60 * 10

logger = logging.getLogger('biostar')

ZERO_COUNTS = dict(
    post_count=0, book_count=0, vote_count=0, mesg_count=0, badge_count=0,
)

SESSION_COUNT_KEY = "counts"

def get_counts(request):
    """
    Returns a populated count object.
    """
    global ZERO_COUNTS

    user = request.user

    if user.is_anonymous():
        # Anonymous users get empty counts.
        return ZERO_COUNTS

    # Authenticated users at this point.
    profile = user.profile

    # No counts found in session.
    missing = SESSION_COUNT_KEY not in request.session

    # How much time has passed since last login.
    elapsed = (utils.now() - profile.last_login).seconds

    # Compute the new counts if necessary.
    if missing or elapsed > settings.SESSION_UPDATE_SECONDS:
        data = dict(
            mesg_count=Message.objects.filter(user=user, unread=True).count(),
            vote_count=Vote.objects.filter(post__author=user,
                                           type__in=(Vote.BOOKMARK, Vote.UP)).count(),
            new_vote_count=Vote.objects.filter(post__author=user,
                                               type__in=(Vote.BOOKMARK, Vote.UP),
                                               unread=True).exclude(author=user).count(),
            post_count=Post.objects.filter(author=user).count(),
            book_count=Vote.objects.filter(author=user, type=Vote.BOOKMARK).count(),
            award_count=Award.objects.filter(user=user).count(),
        )

        # Store the new counts.
        request.session[SESSION_COUNT_KEY] = data

        # Update last login date.
        profile.last_login = utils.now()
        profile.save()

    # Obtain the latest counts.
    counts = request.session.get(SESSION_COUNT_KEY) or ZERO_COUNTS

    return counts


def extras(request):
    # These values will be added to each context

    context = {
        "BIOSTAR_VERSION": VERSION,
        "user": request.user,
        "request": request,
        "recaptcha": settings.RECAPTCHA_PUBLIC_KEY,
        "counts": get_counts(request),
        "TRAFFIC": get_traffic(request),
        "recent_votes": query.recent_votes(request),
        "recent_users": query.recent_users(request),
        "recent_awards": query.recent_awards(request),
        "recent_replies": query.recent_replies(request),
        "SITE_LOGO": settings.SITE_LOGO,
        "SUBDOMAIN_FLAG": request.site.id != settings.SITE_ID
    }

    return context
