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


def add_user_attributes(user):
    """
    Mutates the user in the request  to fill in required attributes.
    """
    if not user.is_authenticated():
        user.is_moderator = user.is_admin = False
    else:
        user.is_admin = (user.type == User.ADMIN)
        user.is_moderator = user.is_admin or (user.type == User.MODERATOR)


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
    parts = text.split(",")
    parts = map(strip, parts)
    parts = filter(None, parts)
    return list(parts)


@transaction.atomic
def create_toplevel_post(data, user, file=None):
    # Creating a top level post from  data
    title = data.get('title', '').strip()
    type = data.get('type', '') or Post.QUESTION
    tag_val = data.get('tag_val', '')
    site = data.get('site', settings.SITE_ID)
    content = data.get('content', '')

    tags = tag_split(tag_val)

    # Add an extra tag if no a question.
    if type != Post.QUESTION:
        tags.append(Post.TYPE_CHOICES_MAP.get(type, ''))

    tag_val = ",".join(tags)

    # Create the post.
    post = Post.objects.create(content=content, title=title, tag_val=tag_val,
                               author=user, type=type, file=file, site_id=site)


    # Return the updated object, otherwise the foreign keys are unset.
    post = Post.objects.get(pk=post.id)

    # Set the tags on the post
    post.tags.set(*tags)

    return post


def can_moderate_post(request, user, post):
    return user.can_moderate_post(post)


def can_moderate_user(request, user, target):
    if target.is_staff:
        return False

    if target.is_admin:
        return False

    if user.is_staff or user.is_moderator:
        return True

    return False


@transaction.atomic
def create_content_post(content, parent, user, post_type=None, file=None):
    # Creating a content level post from data
    post = Post.objects.create(parent=parent, content=content, type=post_type, author=user, file=file)
    return post


def write_access_post(user, post):
    """
    A user may write the post if the post is readable and
    the user is the author of the post or is a moderator
    """
    cond = (user == post.author) or user.is_moderator
    return cond


def write_access_func(user):
    """
    Returns a function that can check a post for
    write access.
    """

    def func(post):
        return write_access_post(user, post)

    return func


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

    return func(request, pk=None, parent=parent)


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