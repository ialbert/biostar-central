import logging
from socket import gethostbyaddr, gethostbyname
from django.conf import settings
from django.core.cache import cache
from django.shortcuts import redirect
from ratelimit.utils import is_ratelimited
from biostar.utils import helpers
from . import util, auth

logger = logging.getLogger("engine")


def limiter(get_response):
    """
    Rate Limiter used to deter anon users
    """

    def middleware(request):
        user = request.user

        # Only check anonymous users
        if user.is_anonymous:
            ip = helpers.get_ip(request)
            triplet = helpers.ip_triplet(request)

            if triplet in settings.WHITELIST_IP:
                return get_response(request)

            # Check if the user should be rate limited within a given time period.
            limited = is_ratelimited(request=request,
                                     fn=get_response,
                                     key=settings.RATELIMIT_KEY,
                                     rate=settings.RATELIMIT_RATE,
                                     increment=True)
            # Redirect to static page if limit reached
            # Might spam the logger view
            # auth.db_logger(text=f'user banned ip={ip}', ipaddr=ip)
            if limited:
                return redirect('/static/message.txt')

        return get_response(request)

    return middleware
