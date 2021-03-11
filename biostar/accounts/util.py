import uuid
from datetime import datetime
from django.utils import timezone
from django.utils.timezone import utc
from django.conf import settings


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def now():
    return timezone.now()
