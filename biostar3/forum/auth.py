"""
Access authorizations are performed here
"""
import os
from string import strip
from .models import *
from django.utils.timezone import utc
from datetime import datetime, timedelta
from django.shortcuts import redirect
from django.contrib import messages
from functools import wraps
from decorator import decorator
from django.contrib.staticfiles import finders


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def ago(hours=0, minutes=0, days=0):
    since = now() - timedelta(days=days, hours=hours, minutes=minutes)
    return since


def tag_split(text):
    lower = lambda x: x.lower() if len(x) > 1 else x
    parts = text.split(",")
    parts = map(strip, parts)
    parts = filter(None, parts)
    parts = map(lower, parts)
    return parts


def create_toplevel_post(data, user, group):
    # Creating a top level post from  data
    title = data.get('title', '').strip()
    type = data.get('type', '')
    tags = data.get('tags', '')
    tags = tag_split(tags)
    content = data.get('content', '')

    # Create the post.
    post = Post.objects.create(content=content, title=title,
                               author=user, type=type, usergroup=group)
    # Set the tags on the post
    post.tags.set(*tags)

    # Self referential ForeignKeys need to be updated explicitly!
    Post.objects.filter(pk=post.pk).update(root_id=post.id, parent_id=post.id)

    # Return the updated object, otherwise the foreign keys are unset.
    post = Post.objects.get(pk=post.id)

    return post


def can_moderate_post(user, post):
    return True

def can_moderate_user(user, target):
    return True

def postsub_get_or_create(user, post, sub_type):
    """
    Gets or creates a postsub for the user
    """
    select, create = PostSub.objects.filter, PostSub.objects.create
    return select(user=user, post=post).first() or create(user=user, post=post, type=sub_type)


def groupsub_get_or_create(user, usergroup, sub_type=None):
    """
    Adds a group sub if it does not exist already.
    """

    # Shortcuts.
    select, create = GroupSub.objects.filter, GroupSub.objects.create

    # Remove the group subscription.
    if sub_type == settings.LEAVE_GROUP:
        return select(user=user, usergroup=usergroup).delete()

    # Is there any subscription for the user and group.
    exists = select(user=user, usergroup=usergroup).first()

    # If there is a subscription and any type will do.
    if exists and not sub_type:
        return exists

    # If exists and it is of the requested type.
    if exists and exists.type == sub_type:
        return exists

    if exists:
        # There is a subscription but needs changing.
        exists.type = sub_type
        exists.save()
        return exists
    else:
        # Create a new subscription
        newsub = create(user=user, usergroup=usergroup, type=settings.DEFAULT_MESSAGES)
        return newsub


def create_content_post(content, parent, post_type, user):
    # Creating a content level post from data
    post = Post.objects.create(parent=parent, content=content, type=post_type, author=user)
    return post


def read_access_post(user, post):
    """
    A user may read the post if the post is in a public group or
    the user is part of the group that the post was made in.
    """
    return post.root.usergroup.public or user.groupsubs.filter(id=post.userroot.group.id).exists()


def write_access_post(user, post):
    """
    A user may write the post if the post is readable and
    the user is the author of the post or is a moderator
    """
    write_cond = (user == post.author) or user.is_moderator
    return write_cond and read_access_post(user=user, post=post)


def thread_write_access(user, root):
    """
    Thread write access check.
    """
    read_cond = read_access_post(post=root, user=user)

    def validator(user, post):
        write_cond = (user == post.author) or user.is_moderator
        return read_cond and write_cond

    return validator


@decorator
def valid_user(func, request, pk, target=None):
    """
    Valid user check.
    """

    if int(pk) == 0:
        return redirect(reverse("account_login"))

    target = User.objects.filter(pk=pk).first()

    if not target:
        messages.error(request, "User with id=%s does not exist." % pk)
        return redirect(reverse("home"))

    return func(request=request, pk=None, target=target)


@decorator
def post_view(func, request, pk, post=None, user=None):
    """
    Post read check.
    """

    user = request.user
    post = Post.objects.filter(pk=pk).first()
    home = redirect(reverse("home"))

    if not post:
        # Post does not exists.
        messages.error(request, "Post with id=%s not found." % pk)
        return home

    if not read_access_post(user=user, post=post):
        # Post exists but may not be read by the user.
        messages.error(request, "This post my not be accessed by this user.")
        return home

    return func(request, pk=None, post=post, user=user)


@decorator
def post_edit(func, request, pk, post=None, user=None):
    """
    Post edit check.
    """

    user = request.user
    post = Post.objects.filter(pk=pk).first()
    home = redirect(reverse("home"))

    if not post:
        # Post does not exists.
        messages.error(request, "Post with id=%s not found" % pk)
        return home

    if not write_access_post(user, post):
        # Post exists but it is not writeable by the user.
        messages.error(request, "Post may not be edited by this user!")
        return home

    return func(request, pk=None, post=post, user=user)


@decorator
def content_create(func, request, pk, parent=None):
    """
    Content create check.
    """

    user = request.user
    error = redirect(reverse("home"))

    parent = Post.objects.filter(pk=pk).select_related("group", "group__groupinfo").first()

    if not parent:
        messages.error(request, "Parent post does not exist. Perhaps it has been deleted.")
        return error

    if not read_access_post(user=user, post=parent):
        messages.error(request, "The thread may not not be accessed by this user!")
        return error

    return func(request, pk=None, parent=parent)


@decorator
def group_access(func, request, pk, group=None, user=None):
    """
    Group access check.
    """

    user = request.user
    group = UserGroup.objects.filter(pk=pk).first()
    error = redirect(reverse("home"))

    if not group:
        # Group does not exists.
        messages.error(request, "Group with id=%s does not exist." % pk)
        return error

    return func(request, pk, group=group, user=user)


@decorator
def group_edit(func, request, pk, group=None, user=None):
    """
    Group edit check.
    """

    user = request.user
    group = UserGroup.objects.filter(pk=pk).first()
    error = redirect(reverse("group_list"))

    if not group:
        # Group does not exists.
        messages.error(request, "Group with id=%s does not exist." % pk)
        return error

    perm = GroupPerm.objects.filter(user=user, usergroup=group, role=GroupPerm.ADMIN).first()

    if not perm:
        # Admin users may edit a group.
        messages.error(request, "Only admin users may edit the group")
        return error

    return func(request=request, pk=None, group=group, user=user)


@decorator
def group_create(func, request, user=None):
    """
    Group create check.
    """

    user = request.user
    error = redirect(reverse("group_list"))

    # How many groups has the user created.
    group_count = UserGroup.objects.filter(owner=user).count()

    if user.score < settings.GROUP_MIN_SCORE:
        # Only users above a treshold in score may create a group.
        messages.error(request, "You need %s reputation to create a group." % (settings.GROUP_MIN_SCORE * 10))
        return error

    if group_count >= settings.GROUP_COUNT_PER_USER:
        # Too many groups created by this user.
        messages.error(request, "Only %s groups may be created by each user." % settings.GROUP_COUNT_PER_USER)
        return error

    return func(request=request, user=user)


def remote_ip(request, key='REMOTE_ADDR'):
    """
    Retrieve the IP address from the request. Does not validate the result.
    """
    # Frontend must set the header correctly.
    ip = request.META.get(key, '0.0.0.0')

    # ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    # ip = ip1 or ip2 or '0.0.0.0'

    return ip


def safe_int(text, default=0, maxval=10000):
    try:
        value = int(text)
        return min(value, maxval)
    except ValueError:
        return default

def safe_remove(path):
    path = os.path.abspath(path)
    path = finders.find(path)
    if path and os.path.isfile(path):
        os.remove(path)