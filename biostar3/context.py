from __future__ import absolute_import, division, print_function, unicode_literals

import logging

from django.conf import settings
from django.core.cache import cache
from biostar3 import VERSION
from django.core.cache import caches
from biostar3.forum.models import Post, Vote, Message
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

def reset_cache(request, key):
    counts = get_counts(request)
    counts[key] = 0
    get_counts(request, counts=counts)


def get_counts(request, counts=None):
    """
    Returns a dictionary with
    """
    if request.user.is_anonymous():
        return {}

    user = request.user
    count_key = COUNT_KEY_PATT % user.id

    if counts:
        # Allows updating the cache as well.
        # This will restart the expiration timer so make it shorter.
        cache.set(count_key, counts, USER_SESSION_TIMEOUT)

    counts = cache.get(count_key)

    if not counts:
        # Counts not found need to be recreated.
        logger.info("hitting the cache %s" % count_key)

        # Save the last login time. Counts are computed relative to that.
        last_login = user.profile.last_login

        # Update the last login field.
        user.profile.last_login = now()
        user.profile.save()

        post_count = Post.objects.filter(author=user).count()
        book_count = Vote.objects.filter(author=user, type=Vote.BOOKMARK).count()
        vote_count = Vote.objects.filter(post__author=user, type__in=(Vote.BOOKMARK, Vote.UP), date__gt=last_login).count()
        mesg_count = Message.objects.filter(user=user, unread=True).count()

        counts = dict(
            post_count=post_count,
            book_count=book_count,
            vote_count=vote_count,
            mesg_count=mesg_count,
            badge_count=0,
        )
        cache.set(count_key, counts, USER_SESSION_TIMEOUT)

    return counts


def shortcuts(request):
    # These values will be added to each context

    context = {
        "BIOSTAR_VERSION": VERSION,
        "user": request.user,
        "group": request.group,
        "request": request,
        "recaptcha": settings.RECAPTCHA_PUBLIC_KEY,
        "counts": get_counts(request),
    }

    return context


