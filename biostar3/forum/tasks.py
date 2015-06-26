from __future__ import absolute_import, division, print_function, unicode_literals

import json
from itertools import *

from celery import shared_task

from .models import *
from . import auth, mailer
from biostar3.utils.compat import *


@shared_task
def compute_user_flair(user):
    """
    This is a fairly compute intensive task. Also the
    flair does not change noticably in short periods of time.
    """


@shared_task
def add_user_location(ip, user):
    """
    Attempts to fill in missing user location.
    """

    if not user.profile.location:
        try:
            url = "http://api.hostip.info/get_json.php?ip=%s" % ip
            logger.info("%s, %s, %s" % (ip, user, url))
            fp = urlopen(url, timeout=3)
            data = fp.read().decode("utf-8")
            data = json.loads(data)
            fp.close()
            location = data.get('country_name', '').title()
            if "unknown" not in location.lower():
                user.profile.location = location
                user.profile.save()
        except KeyError as exc:
            logger.error(exc)


@shared_task
def notify_user(user, post):
    pass

@shared_task
def update_post_subscriptions(post):
    """
    Creates post subscriptions for everyone that should have one
    """
    root = post.root

    if post.is_toplevel:
        # Post tags only exists at top level.
        # Subscribe everyone that watches the tag and does not already have a subscription.
        tag_names = auth.tag_split(post.tag_val)
        followers = User.objects.filter(profile__tags__name__in=tag_names, postsub__isnull=True).distinct()
        PostSub.bulk_insert(post=root, users=followers)

    # Find everyone that has been tagged in the body of the post
    # and does not already have a subscription.
    tagged_names = html.find_tagged_names(post)
    if tagged_names:
        tagged_users = User.objects.filter(profile__tags__name__in=tagged_names, postsub__isnull=True).distinct()
        PostSub.bulk_insert(post=root, users=tagged_users)

    # Manage author subscription.
    exists = PostSub.objects.filter(post=root, user=post.author)

    if not exists:
        prefs = post.author.profile.message_prefs
        if prefs == settings.SMART_MODE:
            if post.is_toplevel:
                prefs = settings.EMAIL_TRACKER
            else:
                prefs = settings.LOCAL_TRACKER

        if prefs in (settings.EMAIL_TRACKER, settings.LOCAL_TRACKER):
            PostSub.objects.create(post=root, user=post.author, type=prefs)

@shared_task
def send_notifications(post):
    """
    Generates messages and emails on a post.
    Invoked on each post creation. Will send emails in an internal loop.
    """

    def wants_email(user):
        # Selects users that have chosen email notifications.
        return not (user.profile.message_prefs in settings.LOCAL_TRACKER)

    # Find all users that have any subscription to a post.
    local_tracker = User.objects.filter(postsub__post=post.root).exclude(postsub__type=settings.NO_MESSAGES).select_related("user")
    local_targets = set(local_tracker)

    # Some of these users will also need to get an email.
    email_tracker = User.objects.filter(postsub__post=post.root, postsub__type=settings.EMAIL_TRACKER).select_related("user")
    email_tracker = chain(email_tracker, User.objects.filter(profile__message_prefs=settings.MAILING_LIST))
    email_targets = set(email_tracker)

    # Never send local messages for the post author.
    if post.author in local_targets:
        local_targets.remove(post.author)

    # Emails should not be  sent to author in SMART_mode.
    if post.author.profile.message_prefs == settings.SMART_MODE:
        if post.author in email_targets:
            email_targets.remove(post.author)

    mailer.post_notifications(post, local_targets=local_targets,
                              email_targets=email_targets)

@shared_task
def add(x, y):
    return x + y