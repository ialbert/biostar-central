from __future__ import absolute_import, division, print_function, unicode_literals

from django.conf import settings
from biostar3.forum import models
from biostar3.forum.models import Post, Vote, Group
from django.contrib.auth import get_user_model
from django.contrib import messages

User = get_user_model()

DEFAULT_GROUP = Group.objects.filter(name=settings.DEFAULT_GROUP_NAME).first()

def recent_votes():
    votes = Vote.objects.filter(post__status=Post.OPEN).select_related("post").order_by("-date")[
            :settings.RECENT_VOTE_COUNT]
    return votes


def get_recent_users():
    users = User.objects.all().select_related("profile").order_by("-profile__last_login")[:settings.RECENT_USER_COUNT]
    return users


def get_toplevel_posts(user, group):
    "Returns posts"
    posts = Post.objects.filter(type__in=Post.TOP_LEVEL, group=group)

    if not user.is_moderator:
        posts = posts.exclude(status=Post.DELETED)

    posts = posts.select_related("root", "author", "lastedit_user").prefetch_related("tags").defer("content", "html")

    return posts


def group_filter(self, name):
    "Performs a query to return posts that belong to a group"
    posts = Post.objects.filter(type__in=Post.TOP_LEVEL, group__name=name)
    posts = posts.select_related("root", "author", "lastedit_user").prefetch_related("tag_set")
    return posts


def get_thread(root, user):
    # Populate the object to build a tree that contains all posts in the thread.
    posts = Post.objects.filter(root=root)
    if not user.is_moderator:
        posts = posts.exclude(status=Post.DELETED)

    posts = posts.select_related("root", "author", "lastedit_user", "author__profile")
    posts = posts.order_by("type", "-has_accepted", "-vote_count", "creation_date")

    return posts