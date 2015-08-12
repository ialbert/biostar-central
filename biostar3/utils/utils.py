__author__ = 'ialbert'

from django.utils.timezone import utc
from datetime import datetime, timedelta

def now():
    return datetime.utcnow().replace(tzinfo=utc)

def ago(hours=0, minutes=0, days=0):
    since = now() - timedelta(days=days, hours=hours, minutes=minutes)
    return since


