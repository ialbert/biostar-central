from django.db.models.signals import post_save
from biostar.apps.main import signals
from biostar.apps.posts.models import Post

# Connect new post creation to actions
post_save.connect(signals.post_create_subscriptions, sender=Post, dispatch_uid="create_subs")

# Creates a messsage to everyone involved
post_save.connect(signals.post_create_messages, sender=Post, dispatch_uid="create_messages")
