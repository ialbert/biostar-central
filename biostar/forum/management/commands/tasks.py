import os
from datetime import timedelta
import random
from datetime import datetime
from biostar import VERSION
from django.core.management.base import BaseCommand
from biostar.forum import models, util, tasks
from biostar.forum.models import Post
from django.conf import settings
from biostar.accounts.models import User
import logging

logger = logging.getLogger('engine')

BACKUP_DIR = os.path.join(settings.BASE_DIR, 'export', 'backup')

BUMP, UNBUMP, AWARD = 'bump', 'unbump', 'award'
CHOICES = [BUMP, UNBUMP, AWARD]


def bump(uids, **kwargs):
    """
    Set post rank the current timestamp
    """

    top = Post.objects.filter(status=Post.OPEN, is_toplevel=True).order_by("-rank")[:10]

    top = random.choice(list(top))
    rank = top.rank
    user = User.objects.filter(is_superuser=True).first()

    if uids:
        uids = uids.split(',')
        Post.objects.filter(uid__in=uids).update(rank=rank, lastedit_user=user)
        logger.debug(f'uids={uids} bumped')

    else:
        # Pick a week from a relatively long time ago.
        # TODO: favor unanswered posts ( or popular post )
        weeks = random.randint(200, 700)
        trange = util.now() - timedelta(weeks=weeks)

        # Pick a random toplevel post within date range.
        post = Post.objects.filter(lastedit_date__gt=trange, is_toplevel=True).order_by('?').first()
        user = User.objects.filter(is_superuser=True).first()

        Post.objects.filter(uid=post.uid).update(rank=rank, lastedit_user=user)

        logger.debug(f'post {post.uid}="{post.title}" bumped')

    return


def unbump(uids, **kwargs):
    """
    Set post rank to creation date
    """
    uids = uids.split(',')
    posts = models.Post.objects.filter(uid__in=uids)
    for p in posts:
        p.rank = p.creation_date
        p.save()


def awards(limit=50, **kwargs):
    """
    Give user awards using a batch method.
    """

    tasks.batch_create_awards(limit=limit)

    return


class Command(BaseCommand):
    help = 'Preform action on list of posts.'

    def add_arguments(self, parser):
        parser.add_argument('--uids', '-u', type=str, required=False, default='', help='List of uids')
        parser.add_argument('--action', '-a', type=str, required=True, choices=CHOICES, default='',
                            help='Action to take.')
        parser.add_argument('--user', dest='pg_user', default="www", help='postgres user default=%(default)s')
        parser.add_argument('--prog', dest='prog', default="/usr/bin/pg_dump",
                            help='the postgres program default=%(default)s')
        parser.add_argument('--outdir', dest='outdir', default=BACKUP_DIR, help='output directory default=%(default)s')
        parser.add_argument('--hourly', dest='hourly', action='store_true', default=False, help='hourly datadump'),
        parser.add_argument('--limit', dest='limit', type=int, default=100, help='How many users'),

    def handle(self, *args, **options):
        action = options['action']

        opts = {BUMP: bump, UNBUMP: unbump, AWARD: awards}

        func = opts[action]
        # print()
        # models.Award.objects.all().delete()

        func(**options)
