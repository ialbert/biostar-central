from __future__ import absolute_import, division, print_function, unicode_literals

from django.conf import settings
from django.core.cache import cache
from biostar3 import VERSION
from django.core.cache import caches
from biostar3.forum.models import Post, Vote
from django.utils.timezone import utc
from datetime import datetime, timedelta

COUNT_KEY_PATT = "count_key_%s"

# When to refresh user data
USER_SESSION_TIMEOUT = 60 * 15


def now():
    return datetime.utcnow().replace(tzinfo=utc)

def ago(hours=0, minutes=0, days=0):
    since = now() - timedelta(days=days, hours=hours, minutes=minutes)
    return since


def get_counts(request):
    """
    Returns a dictionary with
    """
    if request.user.is_anonymous():
        return {}

    user = request.user
    count_key = COUNT_KEY_PATT % user.id

    counts = cache.get(count_key)
    if not counts:
        last_login = user.profile.last_login
        # Update the last login field.
        user.profile.last_login = now()
        user.profile.save()

        post_count = Post.objects.filter(author=user).count()
        book_count = Vote.objects.filter(author=user, type=Vote.BOOKMARK).count()
        vote_count = Vote.objects.filter(author=user, type__in=(Vote.BOOKMARK, Vote.UP), date__gt=last_login).count()

        counts = dict(
            post_count=post_count,
            book_count=book_count,
            vote_count=vote_count,
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


