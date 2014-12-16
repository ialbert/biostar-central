__author__ = 'ialbert'
from django.conf import settings
from .models import Post, Vote
from django.contrib.auth import get_user_model

User = get_user_model()

def recent_votes():
    votes = Vote.objects.filter(post__status=Post.OPEN).select_related("post").order_by("-date")[:settings.RECENT_VOTE_COUNT]
    return votes

def get_recent_users():
    users = User.objects.all().select_related("profile").order_by("-profile__last_login")[:settings.RECENT_USER_COUNT]
    return users