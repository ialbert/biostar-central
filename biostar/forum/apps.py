import logging

from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings

logger = logging.getLogger('engine')



class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_post, sender=self)
        pass




def init_post(sender,  **kwargs):

    from django.contrib.auth import get_user_model
    from . import auth, models

    User = get_user_model()

    name, email = settings.ADMINS[0]

    user = User.objects.filter(email=email).first()

    # Make a couple of test posts
    test_posts = {"foo":dict(tag_val="test foo", post_type=models.Post.FORUM, content="bar"),
                  "bar":dict(tag_val="test bar", post_type=models.Post.DATA, content="foo")
                  }

    for title, val in test_posts.items():

        post = auth.create_post(title=title, author=user, content=val["content"],
                         tag_val=val["tag_val"],post_type=val["post_type"])

        logger.info(f"Created {title} post of {post.get_type_display()}")

    return