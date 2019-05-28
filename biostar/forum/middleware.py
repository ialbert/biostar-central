import logging
from datetime import timedelta
from django.contrib import messages
from django.contrib.auth import logout
from django.conf import settings

from biostar.accounts.models import Profile
import biostar.accounts.auth as accounts_auth
from biostar.message.models import Message
from .util import now
from . import tasks, auth
from .models import Post, Vote


logger = logging.getLevelName("biostar")


def forum_middleware(get_response):

    def middleware(request):

        user, session = request.user, request.session

        # Anonymous users are left alone
        if user.is_anonymous:
            return get_response(request)

        # Banned and suspended users are not allowed
        if auth.is_suspended(user=user):
            messages.error(request, f"Account is {user.profile.get_state_display()}")
            logout(request)

        # Check the user profile.
        tasks.check_profile.spool(request=request, user=user)

        last_login = user.profile.last_login or user.profile.date_joined
        elapsed = (now() - last_login).total_seconds()

        # Update count information inside session
        if elapsed > settings.SESSION_UPDATE_SECONDS:
            # Set the last login time.
            Profile.objects.filter(user=user).update(last_login=now())

            # Store the counts in the session.
            msgs = Message.objects.filter(recipient=user, unread=True, sent_date__gt=last_login).count()
            votes = Vote.objects.filter(post__author=user, date__gt=last_login).exclude(author=user).count()
            counts = dict(votes=votes, msgs=msgs)
            request.session["counts"] = counts

        response = get_response(request)
        # Can process response here after its been handled by the view

        return response

    return middleware
