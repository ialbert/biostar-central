from __future__ import absolute_import
from celery import shared_task

from .models import *
from . import auth, mailer


@shared_task
def create_messages(post):
    """
    Generates messages and emails on a post. Invoked on post creation.
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

    # The message content on the site is different from the email.
    # Email must contain the correct block
    content = mailer.render_node(template=em.template, data=context, name="message")

    body = MessageBody.objects.create(
        author=post.author, subject=em.subj, content=content, html=em.html,
    )

    # Shortcut
    select = PostSub.objects.filter

    # All subscribers get local messages.
    def message_generator(mb):
        # Authors will not get a message when they post.
        subs = select(post=root).exclude(user=post.author).select_related("user").all()
        date = now()
        for sub in subs:
            yield Message(user=sub.user, body=mb, date=date)

    # Bulk insert for all messages.
    Message.objects.bulk_create(message_generator(body), batch_size=100)

    #
    # Generate the email messages.
    # DEFAULT_MESSAGES will only send email to root author on other people's contributions.
    #
    def token_generator(obj):
        date = now()

        # Find everyone that could get an email.
        # This could (probably) be done in a query but the logic gets a little complicated.
        subs = select(post=root, type__in=settings.EMAIL_MESSAGE_TYPES).exclude(user=root.author).select_related(
            "user").all()

        # Check if author has default messaging.
        smart_sub = select(post=root, user=root.author, type=settings.DEFAULT_MESSAGES).select_related("user").first()

        # The author has default subscription and is not the author of the current post.
        if smart_sub and post.author != root.author:
            # We need to flatten to be able to append to it.
            subs = list(subs) + [smart_sub]

        # Generate an email for all candidate.
        for sub in subs:
            token = auth.make_uuid(size=8)
            em.send(to=[sub.user.email], token=token)
            yield ReplyToken(user=sub.user, post=obj, token=token, date=date)

    # Insert the reply tokens into the database.
    ReplyToken.objects.bulk_create(token_generator(post), batch_size=100)


@shared_task
def add(x, y):
    return x + y