import logging
import time
from functools import wraps
from socket import gethostbyaddr, gethostbyname
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import logout
from django.core.cache import cache
from django.shortcuts import redirect
from biostar.accounts.models import Profile, Message
from biostar.accounts.tasks import detect_location

from biostar.utils import helpers

from . import auth, tasks, const, util
from .models import Vote
from .util import now

logger = logging.getLogger("engine")


def benchmark(get_response):
    """
    Prints the time needed to perform a request.
    """

    def middleware(request):

        # Start timer.
        start = time.time()

        # Performs the request
        response = get_response(request)

        # Elapsed time.
        delta = int((time.time() - start) * 1000)

        # Generate timing message.
        msg = f'time={delta}ms for path={request.path}'

        if delta > 1000:
            ip = helpers.get_ip(request)
            uid = request.user.profile.uid if request.user.is_authenticated else '0'
            #agent = request.META.get('HTTP_USER_AGENT', None)
            logger.warning(f"SLOW: {msg} IP:{ip} uid:{uid}")
        elif settings.DEBUG:
            logger.info(f'{msg}')

        return response

    return middleware


def update_status(user):

    # Update a new user into trusted after a threshold score is reached.
    if (user.profile.state == Profile.NEW) and (user.profile.score > 50):
        user.profile.state = Profile.TRUSTED
        user.save()
        return True

    return user.profile.trusted


def user_tasks(get_response):
    """
    Tasks run for authenticated users.
    """

    def middleware(request):

        user, session = request.user, request.session

        # Views for anonymous users are not analzed further.
        if user.is_anonymous:
            return get_response(request)

        # Banned and suspended will be logged out.
        if auth.is_suspended(user=user):
            messages.error(request, f"Account is {user.profile.get_state_display()}")
            logout(request)

        update_status(user=user)

        # Find out the time since the last visit.
        elapsed = (now() - user.profile.last_login).total_seconds()

        # Update information since the last visit.
        if elapsed > settings.SESSION_UPDATE_SECONDS:

            # Detect user location if not set in the profile.
            ip = helpers.get_ip(request)

            # Detect user location if not set in the profile.
            if not user.profile.location:
                detect_location.spool(ip=ip, user_id=user.id)

            # Set the last login time.
            Profile.objects.filter(user=user).update(last_login=now())

            # Compute latest counts.
            counts = auth.get_counts(user=user)

            # Set the session.
            request.session[settings.SESSION_COUNT_KEY] = counts

            # Trigger award generation.
            tasks.create_user_awards.spool(user_id=user.id)

        # Can process response here after its been handled by the view
        response = get_response(request)

        return response

    return middleware
