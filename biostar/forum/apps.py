import logging, random

from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings

logger = logging.getLogger('engine')


class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        from . import signals
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

        #logger.info("initializing badge %s" % badge)


def init_post(sender,  **kwargs):
    from biostar.accounts.apps import init_users
    from django.contrib.auth import get_user_model
    from . import auth, tasks
    from .models import Post

    tasks.foo_task(1,2,3, a="A", b="B", c="C")

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

    users = []
    for i in range(5):
        email = f"User{i}@lvh.me"
        user, flag = User.objects.get_or_create(email=email)
        users.append(user)

    # Type, title, content
    initial = [
        (Post.BLOG, "A blog post", "This is a blog post"),
        (Post.TUTORIAL, "A tutorial post", "This is a tutorial post."),
        (Post.FORUM, "A forum post", "This is a forum post"),
        (Post.QUESTION, "A question post", "This is a question post")
    ]

    posts = []
    for ptype, title, content in initial:
        author = random.choice(users)
        post = auth.create_post(title=title, author=author, content=content, post_type=ptype)
        posts.append(post)


    # Generate answers
    for count in range(10):
        author = random.choice(users)
        parent = random.choice(posts)
        content = f"Answer number {count}"
        answer = auth.create_post(post_type=Post.ANSWER, parent=parent, content=content, author=author)

    # Generate comments
    for count in range(10):
        author = random.choice(users)
        parents= list(Post.objects.order_by("-id"))
        parent = random.choice(parents)
        content = f"Comment number {count}"
        comment = auth.create_post(post_type=Post.COMMENT, parent=parent, content=content, author=author)


