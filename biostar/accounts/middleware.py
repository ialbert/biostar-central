import logging
from socket import gethostbyaddr, gethostbyname
from django.conf import settings
from django.core.cache import cache
from django.shortcuts import redirect
from ratelimit.utils import is_ratelimited
from . import util

logger = logging.getLogger("biostar")


def lookup_host(name, ip):
    """
    Given host name, return ip address after tyring a DNS reverse lookup
    """

    # Construct unique ID based on presumed ip.
    key = f'hostname-{ip}'

    # Refresh cache every 30 days.
    ttl = 3600 * 30

    if key in cache:
        value = cache.get(key, '')
    else:
        # Do reverse DNS lookup using sockets
        value = gethostbyname(name)
        # Set cache
        cache.set(key, value, ttl)
        logger.info(f'Set host cache for ip={value} domain={name}')

    return value


def domain_whitelisted(request):
    """
    Return True if HOST in request is whitelisted.
    """

    # Bail if no domains are whitelisted
    if not settings.WHITELIST_DOMAIN:
        return False

    try:

        # Get hostname from request
        host = util.get_hostname(request)

        # Get IP address from request
        ip = util.get_ip(request)

        # Find IP for this host using DNS lookup
        # and validate it matches IP found in request.
        ip_matches = (ip == lookup_host(host, ip))

        # Return list of whitelisted domains that match current host.
        allowed = list(filter(lambda d: host.endswith(d), settings.WHITELIST_DOMAIN))

        # The domain is whitelisted
        whitelisted = allowed and ip_matches

        return whitelisted

    except Exception as exc:
        logger.info(f"Error checking domain :{exc}")
        return False


def limiter(get_response):
    """
    Rate Limiter used to deter anon users
    """

    def middleware(request):
        user = request.user

        # Only check anonymous users
        if user.is_anonymous:
            ip = util.ip_triplet(request)

            if (ip in settings.WHITELIST_IP) or domain_whitelisted(request):
                return get_response(request)

            # Check if the user should be rate limited within a given time period.
            limited = is_ratelimited(request=request,
                                     fn=get_response,
                                     key=settings.RATELIMIT_KEY,
                                     rate=settings.RATELIMIT_RATE,
                                     increment=True)
            # Redirect to static page if limit reached
            if limited:
                return redirect('/static/message.txt')

        return get_response(request)

    return middleware
