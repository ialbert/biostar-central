__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError
import os
from optparse import make_option
import random
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    help = 'Runs quick patches over the data'

    option_list = BaseCommand.option_list + (
        make_option('--users', dest='users', action='store_true', default=False, help='patches_users'),
        make_option('--bump', dest='bump', action='store_true', default=False, help='bumps a random post'),
        make_option('--bump_id', dest='bump_id', type=int, help='bumps a specific post'),
    )

    def handle(self, *args, **options):

        if options['users']:
            patch_users()

        if options['bump']:
            bump()

        pk = options['bump_id']
        if pk:
            bump(pk)

def patch_users():
    from biostar.apps.users.models import User, Profile
    from biostar.const import DEFAULT_MESSAGES

    users = Profile.objects.all()
    users.update(message_prefs=DEFAULT_MESSAGES)

def bump(pk=None):
    from biostar.apps.posts.models import Post
    from biostar.apps.users.models import User
    from biostar.const import now

    if not pk:
        query = Post.objects.filter(type=Post.QUESTION, status=Post.OPEN)

        value = random.random()

        if value > 0.75:
            query = query.filter(reply_count=0)

        query = query.values_list("id")

        ids = [ p[0] for p in query ]

        pk = random.choice(ids)

    community = User.objects.get(pk=1)
    post = Post.objects.get(pk=pk)
    logger.info(post.title)

    if not post.is_toplevel:
        logger.warning("post is not at toplevel")
    
    post.lastedit_date = now()
    post.lastedit_user = community
    post.save()


