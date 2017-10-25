import logging
import shutil
from django.core.management.base import BaseCommand
from engine.models import Data


logger = logging.getLogger('engine')


def copy(fname1, fname2):
    # copy content of file1 to file2

    copy_from = Data.objects.filter(path=fname1).first()
    copy_to = Data.objects.filter(path=fname2).first()

    if copy_from and copy_to:
        # make sure the data is "Ready"
        shutil.copyfile(fname1, fname2)
        logger.info(f"copied contents of {fname1} to {fname2}")

    else:
        logger.error(f"Files {fname1} and {fname2} not in database.")
    pass


def add(fname, project_id):

    pass



def unpack(fname):

    # Check if tar or gz then unpack it on the same dir
    to_unpack = Data.objects.filter(path=fname).first()

    if to_unpack.state == Data.READY:
        # then unpack.
        pass

    # should it update the data path after unpacking ?
    pass


class Command(BaseCommand):

    help = 'Creates a project.'

    def add_arguments(self, parser):

        parser.add_argument('--copy', action='store_true', default=False,
                            help="Copy contents of fname1 in database to fname2 (alfname1so in database).")

        parser.add_argument('--add', action='store_true', default=False,
                            help="Add file to a project.")

        parser.add_argument('--unpack', action='store_true', default=False,
                            help="Unpack a tar.gz or gzip file already in data base")

        parser.add_argument('--pid',
                            help="Project id to add to.")

        parser.add_argument('--fname1', action='store_true',
                            help="File to add or copy from.")

        parser.add_argument('--fname2', action='store_true',
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
            copy(fname1=fname1, fname2=fname2)

        elif add_data and pid:
            add(files, pid)

        elif unpack_data:
            unpack(files)