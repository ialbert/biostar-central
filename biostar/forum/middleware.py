import logging
from datetime import timedelta
from django.contrib import messages
from django.contrib.auth import logout
from django.conf import settings

from biostar.accounts.models import Profile
import biostar.accounts.auth as accounts_auth
from .util import now
from . import tasks
from .models import Post, Vote


logger = logging.getLevelName("biostar")


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

        # Check the user profile.
        if tasks.HAS_UWSGI:
            tasks.async_check_profile(request=request, user_id=user.id)
        else:
            accounts_auth.check_user_profile(request=request, user=user)

        last_login = user.profile.last_login or user.profile.date_joined
        elapsed = (now() - last_login).total_seconds()

        # Update count information inside session
        if elapsed > settings.SESSION_UPDATE_SECONDS:
            # Set the last login time.
            Profile.objects.filter(user=user).update(last_login=now())

            # Store the counts in the session.
            votes = Vote.objects.filter(post__author=user, date__gt=last_login).exclude(author=user).count()
            counts = dict(votes=votes)
            request.session["counts"] = counts

        response = get_response(request)
        # Can process response here after its been handled by the view

        return response

    return middleware
