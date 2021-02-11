import logging
import time

from socket import gethostbyaddr, gethostbyname
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import logout
from django.core.cache import cache
from django.shortcuts import redirect
from biostar.accounts.models import Profile, Message
from biostar.accounts.tasks import detect_location

from . import auth, tasks, const, util
from .models import Vote
from .util import now

logger = logging.getLogger("biostar")


def get_ip(request):
    """
    Attempts to extract the IP number from the HTTP request headers.
    """
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'
    return ip


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
            logger.warning(f"\n***\n*** SLOW: {msg}\n***\a")
        else:
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


def domain_is_whitelisted(ip):
    try:
        host = gethostbyaddr(ip)[0]
        return host.endswith(settings.WHITE_LIST_DOMAIN) and (ip == gethostbyname(host))
    except:
        return False


def ban_ip(get_response):
    """

    """
    def middleware(request):
        user = request.user

        if settings.DEBUG:
            return get_response(request)

        if user.is_anonymous:
            oip = get_ip(request)
            ips = oip.split(".")[:-1]
            ip = ".".join(ips)

            if ip in settings.IP_WHITELIST:
                return get_response(request)

            if ip not in cache:
                cache.set(ip, 0, settings.TIME_PERIOD)

            value = cache.get(ip)
            if value >= settings.MAX_VISITS:
                # Raise redirect exception
                if domain_is_whitelisted(oip):
                    cache.set(ip, 0)
                else:
                    now = util.now()
                    message = f"{now}\tbanned\t{ip}\t{oip}\n"
                    logger.error(message)
                    fp = open(settings.BANNED_IPS, "a")
                    fp.write(message)
                    fp.close()
                    return redirect('/static/message.txt')
            else:
                cache.incr(ip)

        return get_response(request)

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

        update_status(user=user)
        
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
            vote_count = Vote.objects.filter(post__author=user, date__gt=user.profile.last_login).exclude(author=user).count()

            # Store the counts into the session.
            counts = dict(message_count=message_count, vote_count=vote_count)

            # Set the session.
            request.session[const.COUNT_DATA_KEY] = counts

            # Trigger award generation.
            tasks.create_user_awards.spool(user_id=user.id)

        # Can process response here after its been handled by the view
        response = get_response(request)

        return response

    return middleware
