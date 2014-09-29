"""
Dumps a postgresql database into a file
"""

import os
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from datetime import datetime
from optparse import make_option
import biostar

import logging

logger = logging.getLogger("command")


def abspath(*args):
    """Generates absolute paths"""
    return os.path.abspath(os.path.join(*args))


class Command(BaseCommand):
    help = 'Dumps the postgresql database into a file.'

    option_list = BaseCommand.option_list + (
        make_option('--hourly', dest='hourly', action='store_true', default=False, help='hourly datadump'),
        make_option('-u', dest='pg_user', default="www", help='postgres user default=%default'),
        make_option('-p', dest='prog', default="/usr/bin/pg_dump", help='the postgres program default=%default'),
        make_option('-o', dest='outdir', default="~/data/", help='output directory default=%default'),
    )

    def handle(self, *args, **options):
        pg_user = options['pg_user']
        prog = options['prog']
        hourly = options['hourly']
        outdir = options['outdir']
        main(pg_user=pg_user, hourly=hourly, prog=prog, outdir=outdir)


def main(pg_user, hourly, prog, outdir):
    # Get the full path to the directory.
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

    db_file = "%s-%s-%s.sql.gz" % (pg_name, biostar.VERSION, tstamp)
    db_file = abspath(outdir, db_file)

    params = dict(
        pg_user=pg_user,
        pg_name=pg_name,
        version=biostar.VERSION,
        db_file=db_file,
        prog=prog,
    )

    #logger.info("saving %(pg_name)s to %(db_file)s" % params)

    cmd = "%(prog)s -Fp -x -O -b -U %(pg_user)s %(pg_name)s | gzip > %(db_file)s" % params

    # Running the command
    logger.info("%s" % cmd)
    os.system(cmd)


if __name__ == '__main__':
    #generate_sitemap()
    pass