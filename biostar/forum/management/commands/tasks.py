import os
from datetime import timedelta
import random
from itertools import islice, count
from datetime import datetime
from biostar import VERSION
from django.core.management.base import BaseCommand
from biostar.forum import models, util, tasks
from biostar.forum.models import Post, Diff
from django.conf import settings
from biostar.accounts.models import User
from biostar.utils.decorators import timeit
import logging
import time

logger = logging.getLogger('engine')

BACKUP_DIR = os.path.join(settings.BASE_DIR, 'export', 'backup')

CHOICES = ['bump', 'unbump', 'award']
BUMP, UNBUMP, AWARD = CHOICES


def bump(uids, **kwargs):
    """
    Set post rank the current timestamp
    """

    top = Post.objects.filter(status=Post.OPEN, is_toplevel=True).order_by("-rank")[2:6]

    top = random.choice(list(top))
    rank = top.rank
    user = User.objects.filter(is_superuser=True).first()

    if not uids:

        # Pick open questions from given time period.
        since = util.now() - timedelta(weeks=10)
        query = Post.objects.filter(type=Post.QUESTION,
                                    status=Post.OPEN,
                                    reply_count=0,
                                    creation_date__gt=since)

        # Get valid uids
        query = query.values_list("uid", flat=True)
        uids = [p for p in query]

        # Pick a random uid
        uids = [random.choice(uids)]
    else:
        uids = uids.split(',')

    Post.objects.filter(uid__in=uids).update(rank=rank, lastedit_user=user)
    logger.debug(f'uids={uids} bumped')

    return


def unbump(uids, **kwargs):
    """
    Set post rank to creation date
    """
    uids = uids.split(',')
    posts = models.Post.objects.filter(uid__in=uids)
    for p in posts:
        models.Post.objects.filter(uid__in=uids).update(rank=p.creation_date.timestamp())
        logger.debug(f'title={p.title} uid={p.uid} unbumped.')


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
        parser.add_argument('--limit', dest='limit', type=int, default=100,
                            help='Limit how many users/posts to process.'),

    def handle(self, *args, **options):
        action = options['action']

        opts = {BUMP: bump, UNBUMP: unbump, AWARD: awards}

        func = opts[action]
        # print()
        # models.Award.objects.all().delete()

        func(**options)
