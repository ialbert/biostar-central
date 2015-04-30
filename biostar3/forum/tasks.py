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
def create_messages(post):
    """
    Generates messages and emails on a post.
    Invoked on each post creation.
    """

    # Subscriptions are relative to the root post.
    root = post.root

    # Add a message body for the new post.
    site = Site.objects.get_current()

    # Full url to the post
    post_url = "%s://%s%s" % (settings.SITE_SCHEME, site.domain, reverse("post_view", kwargs=dict(pk=post.id)))

    # Full url to the user
    user_url = "%s://%s%s" % (settings.SITE_SCHEME, site.domain, reverse("user_view", kwargs=dict(pk=post.author.id)))

    # The context that will be passed to the post create template.
    context = dict(post=post, site=site, scheme=settings.SITE_SCHEME,
                   post_url=post_url, user_url=user_url,
                   slug=post.root.usergroup.domain)

    # This is the body of the message that gets created.
    em = mailer.EmailTemplate("post_created_message.html", data=context)

    # Shortcut to filter post subscriptions
    def select(**kwargs):
        return PostSub.objects.filter(post=root, **kwargs).select_related("user")

    # Users with subscription to the post other than post authors.
    targets = select().exclude(user=post.author)

    # Create local messages to all targets
    em.create_messages(author=post.author, targets=targets)

    #
    # Generate the email messages.
    #
    # The DEFAULT_MESSAGES setting will only send email to root author.
    #
    def token_generator(obj):
        now = right_now()
        # Find everyone that could get an email.
        # This could (probably) be done in a query but the logic gets a little complicated.
        subs = select(type__in=settings.EMAIL_MESSAGE_TYPES).exclude(user=root.author)

        # Check if author has default messaging.
        default_sub = select(user=root.author, type=settings.DEFAULT_MESSAGES)\
            .exclude(user=post.author)

        # Chain the results into one.
        subs = chain(subs, default_sub)

        # Generate an email for all candidate.
        for sub in subs:
            token = auth.make_uuid(size=8)
            em.send_email(to=[sub.user.email], token=token)
            yield ReplyToken(user=sub.user, post=obj, token=token, date=now)

    # Insert the reply tokens into the database.
    ReplyToken.objects.bulk_create(token_generator(post), batch_size=100)


@shared_task
def add(x, y):
    return x + y