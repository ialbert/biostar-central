import uuid
from datetime import datetime
from django.utils import timezone
from django.utils.timezone import utc

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

def now():
    return timezone.now()

def get_ip(request):
    """
    Attempts to extract the IP number from the HTTP request headers.
    """
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'
    return ip