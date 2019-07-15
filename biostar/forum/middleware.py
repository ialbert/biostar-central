import logging
from datetime import timedelta
from django.contrib import messages
from django.contrib.auth import logout
from django.conf import settings

from biostar.accounts.models import Profile, Message
from .util import now

from . import auth
from biostar.accounts.tasks import detect_location
from . import tasks

from .models import Post, Vote


logger = logging.getLevelName("biostar")


def get_ip(request):
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'
    return ip


def forum_middleware(get_response):

    def middleware(request):

        user, session = request.user, request.session

        # Anonymous users are not processed.
        if user.is_anonymous:
            return get_response(request)

        # Banned and suspended users are not allowed
        if auth.is_suspended(user=user):
            messages.error(request, f"Account is {user.profile.get_state_display()}")
            logout(request)

        if (user.profile.state == Profile.NEW) and (user.profile.score > 10):
            user.profile.state = Profile.TRUSTED
            user.save()

        ip = get_ip(request)

        # Detect user location if not set in the profile.
        detect_location.spool(ip=ip, user_id=user.id)

        elapsed = (now() - user.profile.last_login).total_seconds()

        # Update count information inside session
        if elapsed > settings.SESSION_UPDATE_SECONDS:

            # Set the last login time.
            Profile.objects.filter(user=user).update(last_login=now())

            # Store the counts in the session.
            message_count = Message.objects.filter(recipient=user, unread=True).count()

            vote_count = Vote.objects.filter(post__author=user, date__gt=user.profile.last_login).exclude(author=user).count()

            # Save the counts into the session.
            counts = dict(message_count=message_count, vote_count=vote_count)
            request.session["counts"] = counts

        response = get_response(request)

        tasks.create_user_awards.spool(user_id=user.id)
        # Can process response here after its been handled by the view

        return response

    return middleware
