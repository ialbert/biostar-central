"""

There are no database models declarations in this file. Data models are specified in the apps.

Only signals and connections between models are specfied here.
"""
from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django.db.models import signals
from biostar.apps.posts.models import Post
from biostar.apps.messages.models import Message, MessageBody
from biostar.apps.posts.models import Subscription

logger = logging.getLogger(__name__)

NEW_POST_CREATED_MESSAGE_TEMPLATE = "messages/post.created.html"

def post_create_subscriptions(sender, instance, created, *args, **kwargs):
    "The actions to undertake when creating a new post"
    from biostar.apps.posts.models import Subscription

    if created:
        # Create a subscription by the user
        Subscription.create(post=instance, user=instance.author)

def post_create_messages(sender, instance, created, *args, **kwargs):
    "The actions to undertake when creating a new post"

    if created:
        # We need to load the module here to avoid circular imports
        # withouth having to introduce a new module.
        from biostar.apps.util import html

        # The user sending the notifications.
        author = instance.author

        # Get all subscriptions for the post.
        subs = Subscription.objects.get_subs(instance).exclude(user=author)

        # Generate the message from the template.
        text = html.render(name=NEW_POST_CREATED_MESSAGE_TEMPLATE, post=instance, user=author)

        # Create the message body.
        body = MessageBody(author=author, subject=instance.title, text=text)
        body.save()

        # This generator will produce the messages.
        def messages():
            for sub in subs:
                yield Message(user=sub.user, body=body)

        # Bulk insert of all messages.
        Message.objects.bulk_create(messages(), batch_size=100)

# Connect new post creation to actions
signals.post_save.connect(post_create_subscriptions, sender=Post, dispatch_uid="create_subs")

# Creates a messsage to everyone involved
signals.post_save.connect(post_create_messages, sender=Post, dispatch_uid="create_messages")
