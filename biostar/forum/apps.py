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


def init_messages(sender, **kwargs):

    from django.contrib.auth import get_user_model
    from . import auth, models

    User = get_user_model()

    name, email = settings.ADMINS[0]

    sender = User.objects.filter(email=email).first()
    if not sender:
        sender = User.objects.create(email=email, username="admin", is_superuser=True)
        sender.set_password(email)

    body = "Hello from the biostar-engine developers, we hope you enjoy the website."
    subject = "Welcome to the biostar-engine!"

    test_2 = User.objects.filter(username="tested").first()
    if not test_2:
        # Create user and send message once.
        test_2 = User.objects.create(username="tested", email="tested@tested")
        recipient_list = [sender, test_2]
        msg = auth.create_local_messages(body=body, subject=subject, rec_list=recipient_list,
                                         sender=sender)

        # Test with a message tree whenever debugging
        if settings.DEBUG:
            msg1 = msg[1]
            msg2 = msg[0]
            models.Message.objects.filter(pk=msg2.pk).update(parent_msg=msg1)
