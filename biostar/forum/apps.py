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




def init_messages(sender, **kwargs):



    return








def init_post(sender,  **kwargs):

    from django.contrib.auth import get_user_model
    from . import auth, models

    User = get_user_model()

    name, email = settings.ADMINS[0]

    user = User.objects.filter(email=email).first()

    # Make a couple of test posts
    blog_title = "Welcome to Biostar-Engine!"
    blog_content = "A small description on the biostar-engine and its use"

    tutorial_title = "Get started with your first project"
    tutorial_content = "This is a tutorial on the in's and out's of the engine."

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