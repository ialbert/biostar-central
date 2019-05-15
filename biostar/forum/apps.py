import logging

from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings

logger = logging.getLogger('engine')


class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        # Triggered upon app initialization.
        #post_migrate.connect(init_post, sender=self)
        #post_migrate.connect(init_awards, sender=self)

        pass


def init_awards(sender,  **kwargs):
    "Initializes the badges"
    from biostar.forum.models import Badge
    from biostar.forum.awards import ALL_AWARDS

    for obj in ALL_AWARDS:
        badge = Badge.objects.filter(name=obj.name)

        if badge:
            continue
        badge = Badge.objects.create(name=obj.name)

        # Badge descriptions may change.
        if badge.desc != obj.desc:
            badge.desc = obj.desc
            badge.icon = obj.icon
            badge.type = obj.type
            badge.save()

        logger.info("initializing badge %s" % badge)


def init_post(sender,  **kwargs):

    from django.contrib.auth import get_user_model
    from . import auth, models

    User = get_user_model()

    name, email = settings.ADMINS[0]

    user = User.objects.filter(email=email).first()

    # Make a couple of tested posts
    blog_title = "Welcome to Biostar-Engine!"
    blog_content = "A small description on the biostar-engine and its use"

    tutorial_title = "Get started with your first project"
    tutorial_content = "This is a tested on the in's and out's of the engine."

    test_posts = {
                blog_title: dict(post_type=models.Post.BLOG, content=blog_content),
                tutorial_title:dict(post_type=models.Post.TUTORIAL,
                                    content=tutorial_content),
                  }

    for title, val in test_posts.items():

        if models.Post.objects.filter(title=title).exists():
            continue

        post = auth.create_post(title=title, author=user, content=val["content"],post_type=val["post_type"])

        logger.info(f"Created {title} post of {post.get_type_display()}")

    return