import logging
import time
import os

from pyftpdlib.filesystems import AbstractedFS
from django.contrib.auth.models import AnonymousUser
from biostar.engine.models import Project, Data, Job
from biostar.engine import auth

from .authorizer import perm_map



logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)


class Bunch(object):
    st_size = 0
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)


def split(path):
    """
    Splits using the os.sep.

    """
    path = [x for x in os.path.normpath(path).split(os.sep) if x]
    return path



def fetch(cond, list, idx):
    # Returns None if cond is not met
    # used to avoid Indexing Errors

    return None if not cond else list[idx]


def query_tab(tab, project, name=None):
    "Return actual path for a path name in /data or /results"

    klass_map = {'data': Data, 'results':Job }
    instance = klass_map[tab].objects.filter(deleted=False, project=project, name=name).first()

    path = ''
    if tab == 'data' and instance:
        path = instance.get_data_dir()
    elif tab == 'results' and instance:
        path = instance.path

    return path



def fetch_file_info(instance, basedir='/'):

    rfname = project = ''
    filetype = 'file'

    if isinstance(instance, Project):
        filetype = "dir"
        rfname = instance.get_project_dir()
        project = instance

    if isinstance(instance, Data):
        rfname = instance.get_data_dir()
        project = instance.project
        filetype = 'file' if len(instance.get_files()) <= 1 else "dir"

    if isinstance(instance, Job):
        rfname = instance.path
        project = instance.project
        filetype = 'dir'

    if isinstance(instance, str):

        path_list = split(basedir)

        project_name = fetch(len(path_list) >= 1, path_list, 0)
        project = Project.objects.filter(name=project_name).first()

        tab = fetch(len(path_list) >= 2, path_list, 1)
        name = fetch(len(path_list) >= 3, path_list, 2)
        base = None if not len(path_list) >= 3 else path_list[3:]

        suffix = instance if not base else os.path.join(*base, instance)
        prefix = query_tab(tab=tab, project=project, name=name)

        rfname = os.path.join(prefix, suffix)

        filetype = 'dir' if rfname and os.path.isdir(rfname) else 'file'


    return project, rfname, filetype



class BiostarFileSystem(AbstractedFS):

    def __init__(self, root, cmd_channel, current_user=None):

        """
         - (str) root: the user "real" home directory (e.g. '/home/user')
         - (instance) cmd_channel: the FTPHandler class instance
         - (str) current_user: username of currently logged in user. Used to key user_table.
        """
        # Set initial current working directory.
        # By default initial cwd is set to "/" to emulate a chroot jail.
        # If a different behavior is desired (e.g. initial cwd = root,
        # to reflect the real filesystem) users overriding this class
        # are responsible to set _cwd attribute as necessary.

        self._cwd = root
        self._root = root
        self.cmd_channel = cmd_channel

        logger.info(f"current_user={current_user}")

        # Dict with all users stored by parent class
        self.user_table = self.cmd_channel.authorizer.user_table

        # Logged in user
        self.user = self.user_table.get(current_user) or dict(user=AnonymousUser)

        # Projects the user has access to
        self.project_list = auth.get_project_list(user=self.user["user"])

        super(BiostarFileSystem, self).__init__(root, cmd_channel)


    def validpath(self, path):
        logger.info(f"path={path}")
        path_list = split(path)

        # Check user did not leave /data or /results while in project.
        tab = fetch(len(path_list) >= 2, path_list, 1)
        if tab and (tab not in ('data', 'results')):
            return False

        return True


    def isdir(self, path):
        logger.info(f"path={path}")
        #return True if path in self.project_list else False
        return True


    def listdir(self, path):
        # This is the root as initialized in the base class

        logger.info(f"path={path}, cwd={self._cwd}")

        path_list = split(path)
        project_name = fetch(len(path_list) >= 1, path_list, 0)
        project = Project.objects.filter(name=project_name).first()
        inside_project = project and len(path_list) == 2

        # List all projects when user is at the root
        if path == self._root:
            return self.project_list

        # We can choose to browse /data or /results inside of a project.
        if project and len(path_list) == 1:
            return ['data', 'results']

        # List project data under /data
        if inside_project and path_list[1] == "data":
            return Data.objects.filter(deleted=False, project=project) or []

        # List project results under /results
        if inside_project and path_list[1] == 'results':
            return Job.objects.filter(deleted=False, project=project) or []

        # At this point, len(path_list) > 2 and
        # the tab should be in (data, results, None ).
        # The latter is guaranteed by self.valid_path().
        name = fetch(len(path_list) >= 3, path_list, 2)

        actual_path = query_tab(tab=path_list[1], project=project, name=name)

        # List the files inside of /data and /results tabs
        if actual_path and len(path_list) == 3:
            # Skip the toc when returning dir contents
            return [os.path.basename(p.path) for p in os.scandir(actual_path)]
        try:
            # List sub-dirs found in /data and /results
            suffix = os.path.join(*path_list[3:])
            full_path = os.path.join(actual_path, suffix)
            return [ os.path.basename(item.path) for item in os.scandir(full_path) ]
        except Exception:
            return []
        

    def chdir(self, path):
        """
        Change the current directory.
        """
        # note: process cwd will be reset by the caller

        logger.info(f"path={path}")
        #self._cwd = f"foo({path})"
        self._cwd = self.fs2ftp(path)


    def format_mlsx(self, basedir, listing, perms, facts, ignore_err=True):

        logger.info(f"basedir={basedir} listing={listing} facts={facts} perms={perms}")
        lines = []
        timefunc = time.localtime if self.cmd_channel.use_gmt_times else time.gmtime
        path_list = split(basedir)

        name = fetch(len(path_list) >= 1, path_list, 0)
        project = Project.objects.filter(name=name).first()

        for instance in listing:

            # Handle /results and /data
            if len(path_list) == 1:
                perm, filetype, rfname = 'elr', 'dir', project.get_project_dir()
            else:
                project, rfname, filetype = fetch_file_info(instance, basedir=basedir)
                perm = perm_map(project=project, user=self.user)

            st = self.stat(rfname)
            unique = "%xg%x" % (st.st_dev, st.st_ino)
            modify = time.strftime("%Y%m%d%H%M%S", timefunc(st.st_mtime))

            lines.append(
                f"type={filetype};size={st.st_size};perm={perm};modify={modify};unique={unique}; {instance}")

        line = "\n".join(lines)
        yield line.encode('utf8', self.cmd_channel.unicode_errors)


