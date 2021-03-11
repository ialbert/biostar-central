import uuid
from datetime import datetime
from django.utils import timezone
from django.utils.timezone import utc
from django.conf import settings


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def now():
    return timezone.now()


def get_ip(request):
    """
    Attempts to extract the IP number from the HTTP request headers.
    """
    # lower the
    key = settings.IP_HEADER_KEY
    meta = request.META

    # Lower case versions of keys
    simple_meta = {k.lower(): v for k, v in request.META.items()}

    ip = meta.get(key, simple_meta.get(key, '0.0.0.0'))

    return ip


def ip_triplet(request):
    """
    Attempt to extract first three number from ip adress.
    """
    oip = get_ip(request=request)
    ips = oip.split(".")[:-1]
    ip = ".".join(ips)
    return ip