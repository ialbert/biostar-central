import os
from datetime import datetime
from biostar import VERSION
from django.core.management.base import BaseCommand
from biostar.forum import models, util, tasks
from django.conf import settings
from biostar.utils.helpers import pg_dump
import logging

logger = logging.getLogger('engine')


BACKUP_DIR = os.path.join(settings.BASE_DIR, 'export', 'backup')

BUMP, UNBUMP, AWARD, DUMP, FAKE = 'bump', 'unbump', 'award', 'pg_dump', 'fake'
CHOICES = [BUMP, UNBUMP, AWARD, DUMP, FAKE]


def fake_awards(**kwargs):
    """
    Creates some fake awards
    """
    #models.Award.objects.all().delete()

    badge = models.Badge.objects.get(name="Autobiographer")
    for user in models.User.objects.all()[:10]:
        award = models.Award(badge=badge, user=user)
        award.save()
        print (f"Awarded {award.badge.name}")
    return

def bump(uids,  **kwargs):
    """
    Set post rank the current timestamp
    """
    uids = uids.split(',')
    rank = util.now().timestamp()
    models.Post.objects.filter(uid__in=uids).update(rank=rank)

    return


def unbump(uids,  **kwargs):
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
        parser.add_argument('--action', '-a', type=str, required=True, choices=CHOICES, default='', help='Action to take.')
        parser.add_argument('--user', dest='pg_user', default="www", help='postgres user default=%(default)s')
        parser.add_argument('--prog', dest='prog', default="/usr/bin/pg_dump", help='the postgres program default=%(default)s')
        parser.add_argument('--outdir', dest='outdir', default=BACKUP_DIR, help='output directory default=%(default)s')
        parser.add_argument('--hourly', dest='hourly', action='store_true', default=False, help='hourly datadump'),
        parser.add_argument('--limit', dest='limit', type=int, default=100, help='How many users'),

    def handle(self, *args, **options):
        action = options['action']

        opts = {BUMP: bump, UNBUMP: unbump, AWARD: awards, DUMP: pg_dump, FAKE:fake_awards}

        func = opts[action]
        #print()
        #models.Award.objects.all().delete()

        func(**options)
