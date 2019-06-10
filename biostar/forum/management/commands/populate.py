import logging
import random

from django.core.management.base import BaseCommand
from django.conf import settings
from django.contrib.auth import get_user_model

from biostar.accounts.apps import init_users
from biostar.forum import auth
from biostar.forum.models import Post


logger = logging.getLogger('engine')

NUSERS = 5
NPOSTS = 10

def init_post(nusers=NUSERS, nposts=NPOSTS):

    # Only initialize when debugging
    if not settings.DEBUG:
        return
    User = get_user_model()

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
    for count in range(nposts):
        author = random.choice(users)
        parent = random.choice(posts)
        content = f"Answer number {count}"
        answer = auth.create_post(post_type=Post.ANSWER, parent=parent, content=content, author=author)

    # Generate comments
    for count in range(nposts):
        author = random.choice(users)
        parents= list(Post.objects.order_by("-id"))
        parent = random.choice(parents)
        content = f"Comment number {count}"
        comment = auth.create_post(post_type=Post.COMMENT, parent=parent, content=content, author=author)


class Command(BaseCommand):
    help = 'Initialize the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--n_users', type=int, default=NUSERS, help="Number of random users to initialize.")
        parser.add_argument('--n_posts', type=int, default=NPOSTS,
                            help="Number of random answers/comments to initialize.")

    def handle(self, *args, **options):

        nusers = options["n_users"]
        nposts = options['n_posts']

        init_users()
        init_post(nposts=nposts, nusers=nusers)

        return
