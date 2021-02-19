from django.core.management.base import BaseCommand
from biostar.forum import models, util, tasks
import logging

logger = logging.getLogger('biostar')

BUMP, UNBUMP, AWARD = 'bump', 'unbump', 'award'

CHOICES = [BUMP, UNBUMP, AWARD]


def bump(uids):
    """
    Set post rank the current timestamp
    """
    rank = util.now().timestamp()
    models.Post.objects.filter(uid__in=uids).update(rank=rank)
    return


def unbump(uids):
    """
    Set post rank to creation date
    """

    posts = models.Post.objects.filter(uid__in=uids)
    for p in posts:
        p.rank = p.creation_date
        p.save()


def awards(limit=10, **kwargs):
    """
    Give user awards using a batch method.
    """

    tasks.batch_create_awards(limit=limit)

    return


class Command(BaseCommand):
    help = 'Preform action on list of posts.'

    def add_arguments(self, parser):
        parser.add_argument('--uids', '-u', type=str, required=False, default='', help='List of uids')
        parser.add_argument('--action', '-a', type=str, required=True, choices=CHOICES,
                            default='', help='Action to take.')

    def handle(self, *args, **options):
        uids = options['uids'].split(',')
        action = options['action']

        opts = {BUMP: bump, UNBUMP: unbump, AWARD: awards}

        func = opts[action]

        #models.Award.objects.all().delete()

        func(uids)
