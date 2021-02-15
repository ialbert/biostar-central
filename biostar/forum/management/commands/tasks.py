from django.core.management.base import BaseCommand
from biostar.forum import models, util
import logging

logger = logging.getLogger('biostar')


BUMP, UNBUMP = 'bump', 'unbump'


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


class Command(BaseCommand):
    help = 'Preform action on list of posts.'

    def add_arguments(self, parser):

        parser.add_argument('--uids', '-u', type=str, default='', help='List of uids')
        parser.add_argument('--action', '-a', type=str, default='', help='Action to take on list of uids')

    def handle(self, *args, **options):

        uids = options['uids']
        action = options['action']

        opts = {BUMP: bump, UNBUMP: unbump}

        if action not in opts:
            logger.error('action={action} not found')
            return

        func = opts[action]
        func(uids)
