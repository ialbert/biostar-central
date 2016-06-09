"""

There are no database models declarations in this file. Data models are specified in the apps.

Only signals and connections between models are specfied here.
"""
from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django.db.models import signals, Q
from allauth.account.signals import user_signed_up
from allauth.socialaccount.signals import social_account_added

from biostar.apps.posts.models import Post, Subscription, ReplyToken
from biostar.apps.messages.models import Message, MessageBody
from biostar.apps.badges.models import Award
from biostar.server.orcid import hook_social_account_added

from biostar.apps.util import html, make_uuid

from django.core import mail
from django.conf import settings
from biostar.const import *
from django.contrib.sites.models import Site

logger = logging.getLogger(__name__)

# This will be the message body on the site.
POST_CREATED_TEXT = "messages/post_created.txt"
POST_CREATED_HTML = "messages/post_created.html"
POST_CREATED_SHORT = "messages/post_created_short.html"

AWARD_CREATED_HTML_TEMPLATE = "messages/award_created.html"

# This will be the message body in an email.

def post_create_messages(sender, instance, created, *args, **kwargs):
    "The actions to undertake when creating a new post"
    from biostar.apps.users.models import User

    post = instance
    if created:
        # The user sending the notifications.
        author = instance.author

        # Insert email subscriptions to users that watch these posts
        if post.is_toplevel:
            cond1 = Q(profile__message_prefs=ALL_MESSAGES)
            cond2 = Q(profile__tags__name__in=post.parse_tags())
            cond = cond1 | cond2
            for watcher in User.objects.filter(cond).exclude(id=author.id):
                sub, flag = Subscription.objects.get_or_create(post=post, user=watcher, type=EMAIL_MESSAGE)

        # Get all subscriptions for the post.
        subs = Subscription.objects.get_subs(post).exclude(user=author)

        # Generate the message from the template.
        content = html.render(name=POST_CREATED_SHORT, post=post, user=author)

        # Generate the email message body.
        site = Site.objects.get_current()
        email_text = html.render(name=POST_CREATED_TEXT, post=post, user=author, site=site)

        # Generate the html message
        email_html = html.render(name=POST_CREATED_HTML, post=post, user=author, site=site)

        # Create the message body.
        body = MessageBody.objects.create(author=author, subject=post.root.title,
                                          text=content, sent_at=post.creation_date)

        # Collects the emails for bulk sending.
        emails, tokens = [], []

        # This generator will produce the messages.
        def messages():
            for sub in subs:
                message = Message(user=sub.user, body=body, sent_at=body.sent_at)

                # collect to a bulk email if the subscription is by email:
                if sub.type in (EMAIL_MESSAGE, ALL_MESSAGES):
                    try:
                        token = ReplyToken(user=sub.user, post=post, token=make_uuid(8), date=now())
                        from_email = settings.EMAIL_FROM_PATTERN % (author.name, settings.DEFAULT_FROM_EMAIL)
                        from_email = from_email.encode("utf-8")
                        reply_to = settings.EMAIL_REPLY_PATTERN % token.token
                        subject = settings.EMAIL_REPLY_SUBJECT % body.subject
                        # create the email message
                        email = mail.EmailMultiAlternatives(
                            subject=subject,
                            body=email_text,
                            from_email=from_email,
                            to=[sub.user.email],
                            headers={'Reply-To': reply_to},
                        )
                        email.attach_alternative(email_html, "text/html")
                        emails.append(email)
                        tokens.append(token)
                    except Exception as exc:
                        # This here can crash the post submission hence the catchall
                        logger.error(exc)

                yield message

        # Bulk insert of all messages. Bypasses the Django ORM!
        Message.objects.bulk_create(messages(), batch_size=100)
        ReplyToken.objects.bulk_create(tokens, batch_size=100)

        try:
            # Bulk sending email messages.
            conn = mail.get_connection()
            conn.send_messages(emails)
        except Exception as exc:
            logger.error("email error %s" % exc)


def award_create_messages(sender, instance, created, *args, **kwargs):
    "The actions to undertake when creating a new post"
    award = instance

    if created:
        # The user sending the notifications.
        user = award.user
        # Generate the message from the template.
        content = html.render(name=AWARD_CREATED_HTML_TEMPLATE, award=award, user=user)

        subject = "Congratulations: you won %s" % award.badge.name

        # Create the message body.
        body = MessageBody.objects.create(author=user, subject=subject, text=content)
        message = Message.objects.create(user=user, body=body, sent_at=body.sent_at)

# Creates a message to everyone involved
signals.post_save.connect(post_create_messages, sender=Post, dispatch_uid="post-create-messages")

# Creates a message when an award has been made
signals.post_save.connect(award_create_messages, sender=Award, dispatch_uid="award-create-messages")


def disconnect_all():
    signals.post_save.disconnect(post_create_messages, sender=Post, dispatch_uid="post-create-messages")
    signals.post_save.disconnect(award_create_messages, sender=Award, dispatch_uid="award-create-messages")


# django-allauth sends a signal when a new user is created using a social provider or a new social
# provider is connected to an existing user.
user_signed_up.connect(hook_social_account_added)
social_account_added.connect(hook_social_account_added)
