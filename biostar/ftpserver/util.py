import logging

import os

from biostar.engine.models import Data, Job

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)



def index(list, idx):
    # Returns None if list is shorter than index
    # used to avoid Indexing Errors and redundancy

    return None if not len(list) >= idx + 1 else list[idx]


def query_tab(tab, project, name=None, show_instance=False):
    "Return actual path for a path name in /data or /results tab"

    klass_map = {'data': Data, 'results':Job }
    instance = klass_map[tab].objects.filter(deleted=False, name=name, project__name=project).first()
    if show_instance:
        return instance

    return None if not instance else instance.get_data_dir()


def parse_virtual_path(ftppath):
    "Parse ftp file path into constituting root_project, tab, pk, and tail. "

    assert isinstance(ftppath, str), ftppath

    # The virtual file path is always seperated by "/" regardless of os.
    path_list = [x for x in os.path.normpath(ftppath).split("/") if x]

    # Project user is under
    root_project = index(path_list, 0) or ''

    # Tab picked ( data or results ).
    tab = index(path_list, 1)

    # Name of specific Data or Job instance.
    instance = index(path_list, 2)

    # Rest of the directory tree to build the actual path with.
    tail = [] if not len(path_list) >= 3 else path_list[3:]

    return root_project, tab, instance, tail
