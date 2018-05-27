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

    auth.create_post(title="foo", author=user, content="bar", tag_val=None,
                    post_type=models.Post.FORUM)

    logger.info("Created test forum post")

    return