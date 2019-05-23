import logging, random

from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings

logger = logging.getLogger('engine')


class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_post, sender=self)
        post_migrate.connect(init_awards, sender=self)


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
    from biostar.accounts.apps import init_users
    from django.contrib.auth import get_user_model
    from . import auth
    from .models import Post

    # Only initialize when debugging
    if not settings.DEBUG:
        return

    init_users()

    User = get_user_model()

    name, email = settings.ADMINS[0]
    user = User.objects.filter(email=email).first()

    # Create admin user.
    if not user:
        user = User.objects.create(email=email, username="admin", is_superuser=True, is_staff=True)
        user.set_password(settings.DEFAULT_ADMIN_PASSWORD)
        user.save()

    # Type, title, content
    initial = [
        (Post.BLOG, "A blog post", "This is a blog post"),
        (Post.TUTORIAL, "A tutorial post", "This is a tutorial post."),
        (Post.FORUM, "A forum post", "This is a forum post"),
        (Post.QUESTION, "A question post", "This is a question post")
    ]

    for ptype, title, content in initial:
        if Post.objects.filter(title=title).exists():
            continue
        post = auth.create_post(title=title, author=user, content=content, post_type=ptype)
        logger.info(f"Created {title} post of {post.get_type_display()}")

    # Generate comments
    for count in range(10):
        all = list(Post.objects.order_by("-id")[:4])
        parent = random.choice(all)
        content = f"Comment number {count}"
        comment = auth.create_post(post_type=Post.COMMENT, parent=parent, content=content, author=user)
