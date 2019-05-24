import logging
import json
from urllib.request import urlopen
from collections import defaultdict

from django.contrib import messages
from django.contrib.auth import logout
from django.conf import settings

from biostar.forum.models import Post, Vote
from biostar.message.models import Message
from biostar.forum.util import now
from biostar.forum import tasks
from biostar.accounts.models import Profile

logger = logging.getLevelName("biostar")


def get_ip(request):
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'
    return ip


def check_user_profile(ip, user):
    logger.info("profile check from %s on %s" % (ip, user))
    if not user.profile.location:
        try:
            url = "http://api.hostip.info/get_json.php?ip=%s" % ip
            logger.info("%s, %s, %s" % (ip, user, url))
            stream = urlopen(url, timeout=3)
            data = json.loads(stream.read())
            stream.close()
            location = data.get('country_name', '').title()
            if "unknown" not in location.lower():
                user.profile.location = location
                user.profile.save()
        except Exception as exc:
            logger.error(exc)


def get_counts(request):
    counts = defaultdict(int)

    # Get the user login,
    user = request.user

    # Compute a few more counts for the user.
    if user.is_authenticated():
        since = user.profile.last_login
        # These are the new messages since the last login.
        user.profile.new_messages = Message.objects.filter(user=user, unread=True, sent_at__gt=since).count()
        user.profile.save()
        # These are the new votes since the last login.
        counts['votes'] = Vote.objects.filter(post__author=user, date__gt=since).count()

    return counts


def forum_middleware(get_response):

    def middleware(request):

        user, session = request.user, request.session

        # Anonymous users are left alone
        if user.is_anonymous:
            return get_response(request)

        # Banned and suspended users are not allowed
        if user.is_authenticated and user.profile.state in (Profile.BANNED, Profile.SUSPENDED):
            messages.error(request, f"Account is {user.profile.get_state_display()}")
            logout(request)
        #
        # # Handle tasks async
        # if tasks.HAS_UWSGI and user.is_authenticated:
        #     tasks.create_user_awards(user_id=user.id)
        #
        # # Handle the user counts and last login.
        # elapsed = (now() - user.profile.last_login).seconds
        #
        # if elapsed > settings.SESSION_UPDATE_SECONDS:
        #     # Set the last login time.
        #     Profile.objects.filter(user=user).update(last_login=now())
        #
        #     # Compute the counts.
        #     counts = get_counts(request)
        #
        #     # Store the counts in the session for later use.
        #     session[settings.SESSION_KEY] = counts
        #
        response = get_response(request)
        # Can process response here after its been handled by the view

        return response

    return middleware
