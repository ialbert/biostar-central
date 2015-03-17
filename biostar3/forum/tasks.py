from __future__ import absolute_import
from celery import shared_task

from .models import *
from . import auth, mailer
from itertools import *
from functools import *
import urllib2, json

@shared_task
def add_user_location(request, user):
    """
    Attempts to fill in missing user location.
    """

    if not user.profile.location:
        try:
            ip = auth.remote_ip(request)
            url = "http://api.hostip.info/get_json.php?ip=%s" % ip
            logger.info("%s, %s, %s" % (ip, user, url))
            f = urllib2.urlopen(url, timeout=3)
            data = json.loads(f.read())
            f.close()
            location = data.get('country_name', '').title()
            if "unknown" not in location.lower():
                user.profile.location = location
                user.profile.save()
        except Exception, exc:
            logger.error(exc)

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

    # The template must contain the blocks subject, message, text and html.
    content = mailer.render_node(template=em.template, data=context, name="message")

    # This is the main message body that will be shown for each message.
    body = MessageBody.objects.create(
        author=post.author, subject=em.subj, content=content, html=em.html,
    )

    # Shortcut to filter post subscriptions
    def select(**kwargs):
        return PostSub.objects.filter(post=root, **kwargs).select_related("user")

    now = right_now()

    # All subscribers get local messages.
    def message_generator(mb):
        # Authors will not get a message when they post.
        subs = select().exclude(user=post.author)

        for sub in subs:
            yield Message(user=sub.user, body=mb, date=now)

    # Bulk insert for all messages.
    Message.objects.bulk_create(message_generator(body), batch_size=100)

    #
    # Generate the email messages.
    #
    # The DEFAULT_MESSAGES setting will only send email to root author.
    #
    def token_generator(obj):

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
            em.send(to=[sub.user.email], token=token)
            yield ReplyToken(user=sub.user, post=obj, token=token, date=now)

    # Insert the reply tokens into the database.
    ReplyToken.objects.bulk_create(token_generator(post), batch_size=100)


@shared_task
def add(x, y):
    return x + y