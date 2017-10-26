import logging
import shutil
import os
import mimetypes
import tarfile
from django.core.management.base import BaseCommand
from engine.models import Data
from biostar.tools.const import COLLECTION_TYPES


logger = logging.getLogger('engine')

def join(*args):
    return os.path.abspath(os.path.join(*args))

#TODO Test
def copy(id1, id2=None, fname=None, project=None):

    # Copies 
    # copy content of file1 to file2

    source = Data.objects.filter(id=id1).first()
    target = Data.objects.filter(id=id2).first()

    if source and target:

        # Make sure the data is "Ready"
        assert source.state == Data.READY
        shutil.copyfile(id1, id2)
        logger.info(f"copied contents of {id1} to {id2}")
        # Change data to ready state
        target.set_ready()

    else:
        logger.error(f"Data {id1} and {id2} not in database.")
    pass


def add(fname, project_id):

    pass


def unzip(fname):


    try:
        mimetype, mimecode = mimetypes.guess_type(fname)
        if mimetype == 'application/x-tar' and mimecode == 'gzip':
            basedir = join(fname, "..")
            os.system(f"cd {basedir} && tar -xvf {fname}")
    except Exception as exc:
        logger.error(f"Error with unzipping {exc}")



def batch_unzip(fnames=()):

    for fname in fnames:
        if not fname:
            continue

        data = Data.objects.filter(rootdir=join(fname, "..")).first()

        if data and data.state == Data.READY:
            assert data.data_type in COLLECTION_TYPES

            unzip(fname)
            logger.info(f"Finished unzipping contents of {fname}.")

    return


class Command(BaseCommand):

    help = 'Creates a project.'

    def add_arguments(self, parser):

        parser.add_argument('--copy', action='store_true', default=False,
                            help="Copy contents of fname1 in database to fname2 (alfname1so in database).")

        parser.add_argument('--add', action='store_true', default=False,
                            help="Add file to a project.")

        parser.add_argument('--unpack', action='store_true', default=False,
                            help="Unpack a tar.gz or gzip file already in data base")

        parser.add_argument('--pid',default=0,
                            help="Project id to add to.")

        parser.add_argument('--fname1',
                            help="File to add or copy from.")

        parser.add_argument('--fname2',
                            help="File to add or copt to.")

    def handle(self, *args, **options):

        fname1 = options.get("fname1")
        fname2 = options.get("fname2")
        pid = options.get("pid")
        copy_data = options.get("copy")
        add_data = options.get("add")
        unpack_data = options.get("unpack")

        if not (fname1 or fname2):
            logger.error("--fname1 or --fname2 need to be set.")

        files = (fname1, fname2)

        if copy_data:
            copy(id1=fname1, id2=fname2)

        elif add_data and pid:
            add(files, pid)

        elif unpack_data:
            batch_unzip(files)
