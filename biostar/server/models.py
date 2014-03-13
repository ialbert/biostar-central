"""

There are no database models declarations in this file. Data models are specified in the apps.

Only signals and connections between models are specfied here.
"""
from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django.db.models import signals

from biostar.apps.users.models import User, Profile
from biostar.apps.posts.models import Post, Subscription
from biostar.apps.messages.models import Message, MessageBody

from biostar.apps.util import html

from django.core import mail
from django.conf import settings
from biostar.const import *
from django.contrib.sites.models import Site

logger = logging.getLogger(__name__)

# This will be the message body on the site.
POST_CREATED_HTML_TEMPLATE = "messages/post.created.html"

# This will be the message body in an email.
POST_CREATED_TEXT_TEMPLATE = "messages/post.created.txt"


def post_create_messages(sender, instance, created, *args, **kwargs):
    "The actions to undertake when creating a new post"

    post = instance
    if created:
        # The user sending the notifications.
        author = instance.author

        # Get all subscriptions for the post.
        subs = Subscription.objects.get_subs(post).exclude(user=author)

        # Generate the message from the template.
        content = html.render(name=POST_CREATED_HTML_TEMPLATE, post=post, user=author)

        # Generate the email message body.
        site = Site.objects.get_current()
        email_text = html.render(name=POST_CREATED_TEXT_TEMPLATE, post=post, user=author, site=site)


        # Create the message body.
        body = MessageBody.objects.create(author=author, subject=post.title,
                                          text=content, sent_at=post.creation_date)

        # Collects the emails for bulk sending.
        emails = []

        # This generator will produce the messages.
        def messages():
            for sub in subs:
                message = Message(user=sub.user, body=body, sent_at=body.sent_at)
                # collect to a bulk email if the subscription is by email:
                if sub.type == EMAIL_MESSAGE:
                    emails.append(
                        (body.subject, email_text, settings.DEFAULT_FROM_EMAIL, [sub.user.email])
                    )
                yield message

            # Generate an email to everyone that has a profile with all messages
            users = User.objects.filter(profile__message_prefs=ALL_MESSAGES)
            for user in users:
                emails.append(
                    (body.subject, email_text, settings.DEFAULT_FROM_EMAIL, [user.email])
                )

        # Bulk insert of all messages. Bypasses the Django ORM!
        Message.objects.bulk_create(messages(), batch_size=100)

        try:
            # Bulk sending email messages.
            results = mail.send_mass_mail(emails)
        except Exception, exc:
            logger.error("email error %s" % exc)

# Creates a message to everyone involved
signals.post_save.connect(post_create_messages, sender=Post, dispatch_uid="post-create-messages")

