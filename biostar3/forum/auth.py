"""
Access authorizations are performed here
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import json
from datetime import timedelta

from django.shortcuts import redirect
from django.contrib import messages
from decorator import decorator

from django.contrib.staticfiles import finders

from .models import *
from biostar3.utils.compat import *


def ago(hours=0, minutes=0, days=0):
    since = right_now() - timedelta(days=days, hours=hours, minutes=minutes)
    return since


def get_group_url(group):
    """
    Find the fully qualified url to a group
    """
    site = Site.objects.get(id=settings.SITE_ID)
    if group.domain == settings.DEFAULT_GROUP_DOMAIN:
        netloc = site.domain
    else:
        netloc = site.domain.split(".")
        if settings.SITE_PREPEND_SUBDOMAIN:
            netloc = [group.domain] + netloc
        else:
            netloc[0] = group.domain
        netloc = ".".join(netloc)

    url = "%s://%s" % (settings.SITE_SCHEME, netloc)
    return url


def tag_split(text):
    lower = lambda x: x.lower() if len(x) > 1 else x
    parts = text.split(",")
    parts = map(strip, parts)
    parts = filter(None, parts)
    parts = map(lower, parts)
    return list(parts)


@transaction.atomic
def create_toplevel_post(data, user, group, file=None):
    # Creating a top level post from  data
    title = data.get('title', '').strip()
    type = data.get('type', '') or Post.QUESTION
    tags = data.get('tags', '')
    tags = tag_split(tags)
    content = data.get('content', '')

    # Create the post.
    post = Post.objects.create(content=content, title=title,
                               author=user, type=type, usergroup=group, file=file)
    # Set the tags on the post
    post.tags.set(*tags)

    # Return the updated object, otherwise the foreign keys are unset.
    post = Post.objects.get(pk=post.id)

    return post


def can_moderate_post(request, user, post):
    if post.author.is_staff:
        return False

    if user.is_staff:
        return True

    if user.is_moderator:
        perm = GroupPerm.objects.filter(user=post.author, usergroup=request.group).first()
        return bool(perm)

    return False


def can_moderate_user(request, user, target):
    if target.is_staff:
        return False

    if user.is_staff:
        return True

    if user.is_moderator:
        # User requests moderation and the target is not a moderator for
        # the current group.
        perm = GroupPerm.objects.filter(user=target, usergroup=request.group).first()
        return bool(perm)

    return False


@transaction.atomic
def postsub_get_or_create(user, post, sub_type):
    """
    Gets or creates a postsub for the user
    """
    select, create = PostSub.objects.filter, PostSub.objects.create
    return select(user=user, post=post).first() or create(user=user, post=post, type=sub_type)


@transaction.atomic
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


@transaction.atomic
def create_content_post(content, parent, user, post_type=None, file=None):
    # Creating a content level post from data
    post = Post.objects.create(parent=parent, content=content, type=post_type, author=user, file=file)
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


CAPTCHA_VERIFY_URL = "https://www.google.com/recaptcha/api/siteverify"

def valid_captcha(request):

    if not settings.RECAPTCHA_PUBLIC_KEY:
        # Captcha validation is not set up.
        return True

    recaptcha_response = request.POST.get('g-recaptcha-response')

    if not recaptcha_response:
        messages.error(request, "Please solve the captcha!")
        return False

    data = {
        'secret': settings.RECAPTCHA_SECRET_KEY,
        'response': recaptcha_response,
        'remoteip': remote_ip(request),
    }

    try:
        # Validate the captcha.
        data = urlencode(data)
        conn = Request(CAPTCHA_VERIFY_URL, data)
        response = urlopen(conn).read()
        result = json.loads(response)
        if not result.get('success'):
            # User failed at solving the capthca.
            messages.error(request, "Failed at the captcha authentication. Please try again!")
            return False
    except Exception as exc:
        # This here is triggered on unexpected errors while solving the capthca.
        logger.error(exc)
        messages.error(request, "Unable to complete captcha challenge: %s" % exc)
        return False

    return True


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