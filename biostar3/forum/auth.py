"""
Access authorizations are performed here
"""
from string import strip
from .models import *
from django.utils.timezone import utc
from datetime import datetime

class AccessDenied(BaseException):
    pass

def now():
    return datetime.utcnow().replace(tzinfo=utc)

def tag_split(text):
    lower = lambda x: x.lower() if len(x) > 1 else x
    parts = text.split(",")
    parts = map(strip, parts)
    parts = filter(None, parts)
    parts = map(lower, parts)
    return parts

def read_access_post(user, post):
    """
    A user may read the post if the post is in a public group or
    the user is part of the group that the post was made in.
    """
    return post.root.group.groupinfo.public or user.groups.filter(name=post.root.group.name).exists()