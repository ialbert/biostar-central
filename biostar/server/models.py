"""

There are no database models declarations in this file. Data models are specified in the apps.

Only signals and connections between models are specfied here.
"""
from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django.db.models import signals
from biostar.apps.posts.models import Post, Subscription
from biostar.apps.util import html
from biostar.apps.messages.models import Message, MessageBody
from django.core import mail
from django.conf import settings
from biostar.const import *

logger = logging.getLogger(__name__)

NEW_POST_CREATED_MESSAGE_TEMPLATE = "messages/post.created.html"

def post_create_messages(sender, instance, created, *args, **kwargs):
    "The actions to undertake when creating a new post"

    post = instance
    if created:
        # The user sending the notifications.
        author = instance.author

        # Get all subscriptions for the post.
        subs = Subscription.objects.get_subs(post).exclude(user=author)

        # Generate the message from the template.
        text = html.render(name=NEW_POST_CREATED_MESSAGE_TEMPLATE, post=post, user=author)

        # Create the message body.
        body = MessageBody(author=author, subject=post.title, text=text)
        body.save()

        # Collects the emails for bulk sending.
        emails = []

        # This generator will produce the messages.
        def messages():
            for sub in subs:
                message = Message(user=sub.user, body=body)
                # collect to a bulk email if the subscription is by email:
                if sub.type == EMAIL_MESSAGE:
                    emails.append(
                        (body.subject, body.text, settings.DEFAULT_FROM_EMAIL, [sub.user.email])
                    )
                yield message

        # Bulk insert of all messages.
        Message.objects.bulk_create(messages(), batch_size=100)


        try:
            # Bulk sending email messages.
            results = mail.send_mass_mail(emails)
        except Exception, exc:
            logger.error("email error %s" % exc)

# Creates a message to everyone involved
signals.post_save.connect(post_create_messages, sender=Post, dispatch_uid="create_messages")
