import logging
import shutil
import os
import tarfile
from django.core.management.base import BaseCommand
from biostar.engine.models import Data, Project
from biostar.tools.const import COLLECTION_TYPES
from biostar.engine import auth

logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


def copy(sourceid=None, targetid=None, fname=None, pid=0):

    project = Project.objects.filter(id=pid).first()

    assert sourceid or fname

    if not targetid:
        # Copies a file to a project
        assert fname and project
        auth.create_data(fname = fname, project=project)
        return

    elif sourceid:
        assert targetid

    source = Data.objects.filter(id=sourceid).first()
    target = Data.objects.filter(id=targetid).first()

    if source and target:

        # Make sure the data is "Ready"
        assert source.state == Data.READY

        # Should the two data types be the same?
        shutil.copyfile(source.file.path, target.file.path)
        logger.info(f"copied contents of {source.name} to {target.name}")
        # Change copied data to ready state
        target.set_ready()

    else:
        logger.error(f"Data {sourceid} and/or {targetid} not in database.")
    return


def add(fname, project_id):

    pass


def unzip(targetid):

    data = Data.objects.filter(id=targetid).first()
    if not data:
        logger.error(f"data.id={targetid} does not exist")
        return

    if data.data_type not in COLLECTION_TYPES:
        logger.error(f"{data.name} wrong data type. Bailing on unzip.")
        return

    fname = data.file.path
    basedir = join(fname, "..")
    logger.info(f"Unpacking data={data.name}, data.file={fname}")
    tar = tarfile.open(fname,"r:gz")
    tar.extractall(path=basedir)
    tar.close()
    logger.info(f"Unzipped data.id={targetid}, name={data.name} to {basedir}")


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

        def method(mode):

            assert target or source
            modemap = {"add":add, "unzip":unzip}
            for data_id in (target, source):
                if data_id:
                    modemap[mode](data_id)

        if not (source or target):
            logger.error("--fname1 or --fname2 need to be set.")

        if copy_data:
            if fname and pid:
                copy(fname=fname, pid=pid)
                return
            copy(sourceid=source, targetid=target)

        elif add_data and pid:
            method("add")

        elif unpack_data:
            method("unzip")
