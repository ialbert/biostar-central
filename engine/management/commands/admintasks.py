import logging
import shutil
import os
import mimetypes
import tarfile
from django.core.management.base import BaseCommand
from engine.models import Data, Project
from biostar.tools.const import COLLECTION_TYPES


logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


def copy(sourceid=None, targetid=None, fname=None, pid=0):

    project = Project.objects.filter(id=pid).first()

    assert sourceid or fname

    if not targetid:
        # Copies a file to a project
        assert fname and project
        project.create_data(fname = fname)
        return

    elif sourceid:
        assert targetid

    source = Data.objects.filter(id=sourceid).first()
    target = Data.objects.filter(id=targetid).first()

    if source and target:

        # Make sure the data is "Ready"
        assert source.state == Data.READY
        shutil.copyfile(source.file.path, target.file.path)
        logger.info(f"copied contents of {source.name} to {target.name}")
        # Change copied data to ready state
        target.set_ready()

    else:
        logger.error(f"Data {sourceid} and {targetid} not in database.")
    return


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

        parser.add_argument('--source',
                            help="Data id to copy from.")

        parser.add_argument('--target',
                            help="Data id to copt to.")

        parser.add_argument('--fname',
                            help="File copy to project.")

    def handle(self, *args, **options):

        source = options.get("source")
        target = options.get("target")
        pid = options.get("pid")
        fname = options.get("fname")
        copy_data = options.get("copy")
        add_data = options.get("add")
        unpack_data = options.get("unpack")

        if not (source or target):
            logger.error("--fname1 or --fname2 need to be set.")

        files = (source, target)

        if copy_data:
            if fname and pid:
                copy(fname=fname, pid=pid)
                return
            copy(sourceid=source, targetid=target)

        elif add_data and pid:
            add(files, pid)

        elif unpack_data:
            batch_unzip(files)
