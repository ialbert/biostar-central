"""
Access authorizations are performed here
"""
from string import strip
from .models import *
from django.utils.timezone import utc
from datetime import datetime, timedelta
from django.shortcuts import redirect
from django.contrib import messages
from functools import wraps

class AccessDenied(BaseException):
    pass

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
                              author=user, type=type, group=group)
    # Set the tags on the post
    post.tags.set(*tags)
    
    # Self referential ForeignKeys need to be updated explicitly!
    Post.objects.filter(pk=post.pk).update(root_id=post.id, parent_id=post.id)

    # Return the updated object, othewise the foreign keys are not set.
    post = Post.objects.get(pk=post.id)

    return post

def create_content_post(data, parent, post_type, user):
    # Creating a content level post from data
    content = data.get('content', '')
    post = Post.objects.create(parent=parent, content=content, type=post_type, author=user)
    return post

def read_access_post(user, post):
    """
    A user may read the post if the post is in a public group or
    the user is part of the group that the post was made in.
    """
    return post.root.group.public or user.usergroups.filter(name=post.root.group.name).exists()

def write_access_post(user, post):
    """
    A user may write the post if the post is readable and
    the user is the author of the post or is a moderator
    """
    write_cond = (user == post.author) or user.is_moderator
    return write_cond and read_access_post(user=user, post=post)

def thread_write_access(user, root):

    read_cond = read_access_post(post=root, user=user)

    def validator(user, post):
        write_cond = (user == post.author) or user.is_moderator
        return read_cond and write_cond

    return validator


def valid_user(function=None):
    """
    Decorator for views that checks that the userid in the pk field is a valid
    user for the current state.
    """

    def decorator(request, pk, *args, **kwargs):
        user = User.objects.filter(pk=pk).first()

        if not user:
            messages.error(request, "User with id=%s not found" % pk)
            return redirect(reverse("home"))

        return function(request, pk, user=user)

    return decorator

def post_read(function=None):
    """
    Decorator for views that checks that the post id in the pk field
    is a valid read target for the current user.
    """

    @wraps(function)
    def decorator(request, pk, *args, **kwargs):
        user = request.user
        post = Post.objects.filter(pk=pk).first()
        home = redirect(reverse("home"))

        if not post:
            # Post does not exists.
            messages.error(request, "Post with id=%s not found" % pk)
            return home

        if not read_access_post(user=user, post=post):
            # Post exists but may not be read by the user.
            messages.error(request, "This post my not be accessed by this user.")
            return home

        return function(request, post=post, user=user)

    return decorator

def post_edit(function=None):
    """
    Decorator for views that checks that the post id in the pk field
    is a valid edit target for the current user.
    """

    @wraps(function)
    def decorator(request, pk, *args, **kwargs):
        user = request.user
        post = Post.objects.filter(pk=pk).first()
        back = redirect(post.get_absolute_url())

        if not write_access_post(user, post):
            # Post exists but it is not writeable by the user.
            messages.error(request, "Post may not be edited by this user!")
            return back

        return function(request, post=post, user=user)

    return decorator

def content_create(function=None):
    """
    Decorator to check content creation on a parent post.
    """

    def decorator(request, pk, *args, **kwargs):
        user = request.user
        error = redirect(reverse("home"))

        parent = Post.objects.filter(pk=pk).select_related("group", "group__groupinfo").first()

        if not parent:
            messages.error(request, "Parent post does not exist. Perhaps it has been deleted.")
            return error

        if not read_access_post(user=user, post=parent):
            messages.error(request, "The thread may not not be accessed by this user!")
            return error

        return function(request, parent=parent,  **kwargs)

    return decorator


def content_write(function=None):
    """
    Decorator to check content creation on a parent post.
    """

    def decorator(request, parent,  **kwargs):
        user = request.user
        error = redirect(reverse("home"))

        if not write_access_post(user=user, post=parent):
            messages.error(request, "This post may not be edited by this user!")
            return error

        return function(request, parent=parent, **kwargs)

    return decorator


def group_access(function=None):
    """
    Decorator for views that checks that the group id the pk field is
    a valid user for the current request.
    """

    def decorator(request, pk, *args, **kwargs):
        user = request.user
        group = UserGroup.objects.filter(pk=pk).first()
        error = redirect(reverse("home"))

        if not group:
            # Group does not exists.
            messages.error(request, "Group with id=%s not found" % pk)
            return error

        return function(request, group=group, user=user, **kwargs)

    return decorator

@group_access
def group_edit(function=None):
    """
    Decorator for views that checks that the group id  the pk field is
    a valid user for the current request.
    """

    def decorator(request, pk, *args, **kwargs):
        user = request.user
        group = UserGroup.objects.filter(pk=pk).first()
        error = redirect(reverse("group_list"))

        if group.owner != user:
            # Only group owners may edit a group.
            messages.error(request, "Only the group owner may edit a group")
            return error

        return function(request, group=group, user=user, **kwargs)

    return decorator

def group_create(function=None):
    """
    Decorator for views that checks that the group id  the pk field is
    a valid user for the current state.
    """

    def decorator(request, *args, **kwargs):
        user = request.user
        error = redirect(reverse("group_list"))
        group_count = UserGroup.objects.filter(owner=user).count()


        if user.score < settings.GROUP_MIN_SCORE:
            # Only users above a treshold in score may create a group.
            messages.error(request, "You need %s reputation to create a group" % (settings.GROUP_MIN_SCORE * 10))
            return error

        if group_count >= settings.GROUP_COUNT_PER_USER:
            # Only group owners may edit a group.
            messages.error(request, "Only %s groups may be created by each user" % settings.GROUP_COUNT_PER_USER)
            return error

        return function(request, user=user)

    return decorator


def remote_ip(request, key='REMOTE_ADDR'):
    """
    Retrieve the IP address from the request.
    Does not validate the result.
    """
    # Frontend must set the header correctly.
    ip = request.META.get(key, '0.0.0.0')

    #ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    #ip = ip1 or ip2 or '0.0.0.0'

    return ip