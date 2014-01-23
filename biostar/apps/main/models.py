from django.db.models import signals
from biostar.apps.main import dispatch
from biostar.apps.posts.models import Post

# Connect new post creation to actions
signals.post_save.connect(dispatch.post_create_subscriptions, sender=Post, dispatch_uid="create_subs")

# Creates a messsage to everyone involved
signals.post_save.connect(dispatch.post_create_messages, sender=Post, dispatch_uid="create_messages")

