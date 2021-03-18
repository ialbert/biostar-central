import logging
import random
import os
from django.core.management.base import BaseCommand
from django.conf import settings
from django.contrib.auth import get_user_model
from biostar.accounts.models import Message, MessageBody
from biostar.forum import auth
from biostar.forum.models import Post, Vote


logger = logging.getLogger('engine')

NUSERS = 0
NPOSTS = 0
User = get_user_model()

def init_post(nusers=NUSERS, nposts=NPOSTS):

    name, email = settings.ADMINS[0]
    user = User.objects.filter(email=email).first()

    # Create admin user.
    if not user:
        user = User.objects.create(email=email, username="admin", is_superuser=True, is_staff=True)
        user.set_password(settings.DEFAULT_ADMIN_PASSWORD)
        user.save()

    users = []
    for i in range(nusers):
        email = f"User{i}@lvh.me"
        user, flag = User.objects.get_or_create(email=email)
        User.objects.filter(pk=user.pk).update(username=f'auto-{user.pk}')
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
        post = Post.objects.create(title=title, author=author, content=content, type=ptype)
        posts.append(post)

    # Drop one post from targets.
    targets = posts[1:]

    # Generate answers
    for count in range(nposts):
        author = random.choice(users)
        parent = random.choice(targets)
        content = f"Answer number {count}"
        answer = Post.objects.create(type=Post.ANSWER, parent=parent, content=content, author=author)

    # Generate comments
    for count in range(nposts):
        author = random.choice(users)
        # Exclude the first generated post.
        parents = list(Post.objects.order_by("-id").exclude(pk=posts[0].pk))
        parent = random.choice(parents)
        content = f"Comment number {count}"
        comment = Post.objects.create(type=Post.COMMENT, parent=parent, content=content, author=author)


def init_messages(nmsgs):
    # Fetch target user to receive messages
    target = User.objects.filter(email=settings.ADMIN_EMAIL).first()
    # Any staff user is used as a source to send messages
    source = User.objects.filter(is_staff=True).first()
    if source and target:
        for m in range(nmsgs):
            body = html = f"This is a test message {m}"
            body = MessageBody.objects.create(body=body, html=html)
            Message.objects.create(sender=source, recipient=target, body=body)

    logger.info(f"Finished initializing messages {nmsgs} messages from:{source} to:{target}")
    return


def init_votes(nvotes):
    # Fetch a post by the target user
    target = User.objects.filter(email=settings.ADMIN_EMAIL).first()
    post = Post.objects.filter(author=target).first()
    if not post:
        post = Post.objects.create(title="Question post", author=target,
                                   content="This is a question post", type=Post.QUESTION)

    # Have source user upvote posts by target user.
    source = User.objects.exclude(pk=target.pk).filter(is_staff=True).first()
    for v in range(nvotes):
        Vote.objects.create(author=source, post=post, type=Vote.UP)

    logger.info(f"Finished initializing {nvotes} up votes from:{source} to:{target} post-uid:{post.uid}")
    return


class Command(BaseCommand):
    help = 'Initialize the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--n_users', type=int, default=NUSERS, help="Number of random users to initialize.")
        parser.add_argument('--n_messages', type=int, default=NUSERS, help="Number of messages to initialize.")
        parser.add_argument('--n_votes', type=int, default=NUSERS, help="Number of votes to initialize.")
        parser.add_argument('--demo', action="store_true", default=False, help="Load demo data")
        parser.add_argument('--n_posts', type=int, default=NPOSTS,
                            help="Number of random answers/comments to initialize.")

    def handle(self, *args, **options):

        nusers = options["n_users"]
        nposts = options['n_posts']
        nmsgs = options['n_messages']
        nvotes = options['n_votes']
        demo = options['demo']

        # Set fields needed for quick demo
        if demo:
            nusers = nposts = 10

        # Only initialize when debugging
        if not settings.DEBUG:
            logger.info("Can not initialize when DEBUG=False")
            return

        logger.info("Populating")
        if nusers or nposts:
            init_post(nposts=nposts, nusers=nusers)
        if nmsgs:
            init_messages(nmsgs=nmsgs)
        if nvotes:
            init_votes(nvotes=nvotes)

        return
