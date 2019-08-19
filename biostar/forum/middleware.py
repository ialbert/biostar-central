import logging
import time

from django.conf import settings
from django.contrib import messages
from django.contrib.auth import logout

from biostar.accounts.models import Profile, Message
from biostar.accounts.tasks import detect_location
from . import auth, tasks, const
from .models import Vote
from .util import now

logger = logging.getLogger("biostar")


def get_ip(request):
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'
    return ip


def benchmark(get_response):
    """
    Benchmarking each request.
    """

    def middleware(request):
        start = time.time()

        response = get_response(request)

        delta = (time.time() - start)*1000

        logger.info(f'time={delta:.0f}ms path={request.path}')

        return response

    return middleware


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

        # Update a new user into trusted after 10 votes.
        # TODO: change to a separate function and a different policy.
        if (user.profile.state == Profile.NEW) and (user.profile.score > 10):
            user.profile.state = Profile.TRUSTED
            user.save()

        # Parses the ip of the request.
        ip = get_ip(request)

        # Find out the time since the last visit.
        elapsed = (now() - user.profile.last_login).total_seconds()

        # Update information since the last visit.
        if elapsed > settings.SESSION_UPDATE_SECONDS:
            # Detect user location if not set in the profile.
            detect_location.spool(ip=ip, user_id=user.id)

            # Set the last login time.
            Profile.objects.filter(user=user).update(last_login=now())

            # The number of new messages since last visit.
            message_count = Message.objects.filter(recipient=user, unread=True).count()

            # The number of new votes since last visit.
            vote_count = Vote.objects.filter(post__author=user, date__gt=user.profile.last_login).exclude(
                author=user).count()

            # Store the counts into the session.
            counts = dict(message_count=message_count, vote_count=vote_count)

            # Set the session.
            request.session[const.COUNT_DATA_KEY] = counts

            tasks.create_user_awards.spool(user_id=user.id)
            # Can process response here after its been handled by the view

        response = get_response(request)

        return response

    return middleware
