import os
from datetime import datetime
from biostar import VERSION
from django.core.management.base import BaseCommand
from biostar.forum import models, util, tasks
from django.conf import settings
import logging

logger = logging.getLogger('biostar')

BUMP, UNBUMP, AWARD, DUMP = 'bump', 'unbump', 'award', 'dump'

CHOICES = [BUMP, UNBUMP, AWARD, DUMP]


def pg_dump(prog, pg_user, outdir, hourly=False,  **kwargs):
    # Get the full path to the directory.
    """
    mac - /usr/local/bin/pg_dump
    linux - /usr/bin/pg_dump
    """

    outdir = os.path.expanduser(outdir)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    pg_name = settings.DATABASE_NAME

    # Hourly database dumps have a simpler name so
    # that they overwrite each other.
    if hourly:
        # These names only include the hours.
        tstamp = datetime.now().strftime("hourly-%H")
    else:
        # These names include the date.
        tstamp = datetime.now().strftime("%Y-%m-%d")

    db_file = "%s-%s-%s.sql.gz" % (pg_name, VERSION, tstamp)
    db_file = os.path.abspath(os.path.join(outdir, db_file))

    params = dict(
        pg_user=pg_user,
        pg_name=pg_name,
        version=VERSION,
        db_file=db_file,
        prog=prog,
    )

    # logger.info("saving %(pg_name)s to %(db_file)s" % params)

    cmd = "%(prog)s -x -O -b -U %(pg_user)s %(pg_name)s | gzip > %(db_file)s" % params

    # Running the command
    logger.info("%s" % cmd)
    os.system(cmd)

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
        parser.add_argument('--action', '-a', type=str, required=True, choices=CHOICES, default='', help='Action to take.')
        parser.add_argument('-u', dest='pg_user', default="www", help='postgres user default=%default')
        parser.add_argument('-p', dest='prog', default="/usr/local/bin/pg_dump", help='the postgres program default=%default')
        parser.add_argument('-o', dest='outdir', default="~/data/", help='output directory default=%default')
        parser.add_argument('--hourly', dest='hourly', action='store_true', default=False, help='hourly datadump'),

    def handle(self, *args, **options):
        action = options['action']

        opts = {BUMP: bump, UNBUMP: unbump, AWARD: awards, DUMP: pg_dump}

        func = opts[action]

        #models.Award.objects.all().delete()

        func(**options)
