import logging
from django.conf import settings
from django.shortcuts import redirect
from ratelimit.utils import is_ratelimited
from . import util

logger = logging.getLogger("biostar")


def limiter(get_response):
    """

    """
    def middleware(request):
        user = request.user

        # Rate limiting is disabled in the settings.
        if not settings.ENABLE_RATELIMIT:
            return get_response(request)

        # Only check anonymous users
        if user.is_anonymous:

            ip = util.ip_triplet(request)
            if ip in settings.IP_WHITELIST:
                return get_response(request)

            # Check if the user should be rate limited
            limited = is_ratelimited(request=request,
                                     fn=get_response,
                                     key=settings.RATELIMIT_KEY,
                                     rate=settings.RATELIMIT_RATE,
                                     increment=True)
            # Redirect to static page.
            if limited:
                return redirect('/static/message.txt')

        return get_response(request)

    return middleware