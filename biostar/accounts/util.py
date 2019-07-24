import uuid
from datetime import datetime

from django.utils.timezone import utc

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

def now():
    return datetime.utcnow().replace(tzinfo=utc)

