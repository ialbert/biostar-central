import logging

import os

from biostar.engine.models import Data, Job

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)



def index(list, idx):
    # Returns None if list is shorter than index
    # used to avoid Indexing Errors and redundancy

    return None if not len(list) >= idx + 1 else list[idx]


def split(path):
    """
    Completely splits path using the os.sep.

    """
    path = [x for x in os.path.normpath(path).split(os.sep) if x]
    return path


def query_tab(tab, pk=None, show_instance=False):
    "Return actual path for a path name in /data or /results tab"

    klass_map = {'data': Data, 'results':Job }
    instance = klass_map[tab].objects.filter(deleted=False, pk=pk).first()
    if show_instance:
        return instance

    return None if not instance else instance.get_data_dir()


def filesystem_mapper(queryset=None, tag=""):
    "Takes queryset returns correct ftp path to be parsed"

    return [f"{p.pk}.{tag}{p.name}" for p in queryset.order_by('pk')] or []


def parse_pk(string):
    "Parse a project pk from given project name."

    try:
        pk = string.split(".")[0].strip()
        return int(pk)
    except Exception as exc:
        logger.error(f"{exc}")
        return None


def parse_virtual_path(ftppath):
    "Parse ftp file path into constituting root_project, tab, pk, and tail. "

    assert isinstance(ftppath, str), ftppath

    path_list = split(ftppath)

    # Project user is under
    root_project = index(path_list, 0) or ''

    # Build filesystem using project pk value since names can the same.
    root_project = parse_pk(string=root_project)

    # Tab picked ( data or results ).
    tab = index(path_list, 1)

    # Pk of specific Data or Job instance.
    instance = index(path_list, 2)

    pk = parse_pk(string=instance)

    # Rest of the directory tree to build the actual path with.
    tail = [] if not len(path_list) >= 3 else path_list[3:]

    return root_project, tab, pk, tail