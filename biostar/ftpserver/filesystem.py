import logging
import time
import os
import sys

from pyftpdlib.filesystems import AbstractedFS
from django.contrib.auth.models import AnonymousUser
from biostar.engine.models import Project, Data, Job
from biostar.engine import auth
from biostar.accounts.models import Profile

from .authorizer import perm_map



logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)




def fetch(cond, list, idx):
    # Returns None if cond is not met
    # used to avoid Indexing Errors

    return None if not cond else list[idx]


def split(path):
    """
    Completely splits path using the os.sep.

    """
    path = [x for x in os.path.normpath(path).split(os.sep) if x]
    return path


def query_tab(tab, project, name=None):
    "Return actual path for a path name in /data or /results tab"

    klass_map = {'data': Data, 'results':Job }
    instance = klass_map[tab].objects.filter(deleted=False,
                                             project=project, name=name).first()

    return None if not instance else instance.get_data_dir()


def fetch_file_info(instance, project=None, tab='data', base=[], name=None):
    """
    Fetch the actual file path and filetype of a given instance.
    """

    rfname = None
    data_file = isinstance(instance, Data) and len(instance.get_files()) <= 1
    filetype = 'file' if data_file else "dir"

    # Fetch files when the instance is a database model.
    if isinstance(instance, Data) or isinstance(instance, Job) or isinstance(instance, Project):
        rfname = instance.get_data_dir()

    # Build the real path when the instance is an ftp file path.
    if isinstance(instance, str):

        suffix = instance if not base else os.path.join(*base, instance)
        prefix = query_tab(tab=tab, project=project, name=name)
        rfname = os.path.join(prefix, suffix)
        is_dir = os.path.exists(rfname) and os.path.isdir(rfname)
        filetype = 'dir' if is_dir else 'file'


    return rfname, filetype



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
        """
        Deconstructs the path into :
            - the root_project
            - current tab
            - Data or Job name
            - rest of directory tree
        Validates each part to make sure it exists and is valid.
        Also checks user is not banned or suspended.

        :param path: current path requested by user
        :return: True if path is valid else False
        """

        logger.info(f"path={path}")

        # Check the user still has access.
        self._check_user_status()
        path_list = split(path)

        root_project, tab, name, base = self.extract(path_list=path_list)
        project = Project.objects.filter(name=root_project)
        data = Data.objects.filter(name=name, project=project.first(), deleted=False)
        results = Job.objects.filter(name=name, project=project.first(), deleted=False)

        # Root project must exist
        if root_project and not project.exists():
            return False

        # Check user did not leave /data or /results while in project.
        if tab and (tab not in ('data', 'results')):
            return False

        # Data or Job must be valid.
        if name and not ( data or results ):
            return False

        # The path is valid at this point
        is_valid = True
        if base:
            actual_path = query_tab(tab=tab, project=project.first(), name=name)
            suffix = os.path.join(*base)
            is_valid = os.path.exists( os.path.join(actual_path, suffix))

        logger.info(f"path={is_valid}")
        return is_valid



    def isdir(self, path):
        logger.info(f"path={path}")
        #return True if path in self.project_list else False
        return True


    def listdir(self, path):
        # This is the root as initialized in the base class

        logger.info(f"path={path}, cwd={self._cwd}")

        path_list = split(path)
        root_project, tab, name, base = self.extract(path_list=path_list)

        project = Project.objects.filter(name=root_project).first()
        inside_project = project and len(path_list) == 2

        # List all projects when user is at the root
        if path == self.root:
            #return [f"Project-{p}" for p in self.project_list]
            return self.project_list

        # We can choose to browse /data or /results inside of a project.
        if project and len(path_list) == 1:
            return ['data', 'results']

        # List project data under /data
        if inside_project and tab == "data":
            return Data.objects.filter(deleted=False, project=project) or []

        # List project results under /results
        if inside_project and tab == 'results':
            return Job.objects.filter(deleted=False, project=project) or []

        # At this point, len(path_list) > 2 and
        # the tab should be in (data, results).
        # The latter is guaranteed by self.valid_path.
        actual_path = query_tab(tab=tab, project=project, name=name)

        # List the files inside of /data and /results tabs
        if actual_path and len(path_list) == 3:
            return [os.path.basename(p.path) for p in os.scandir(actual_path)]

        try:
            # List sub-dirs found in /data and /results
            suffix = os.path.join(*base)
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

        root_project, tab, name, base = self.extract(path_list=path_list)

        project = Project.objects.filter(name=root_project).first()

        for instance in listing:

            # Handle /results and /data
            if len(path_list) == 1:
                perm, filetype, rfname = 'elr', 'dir', project.get_data_dir()
            else:
                rfname, filetype = fetch_file_info(instance, tab=tab, name=name,
                                                   project=project, base=base)
                perm = perm_map(project=project, user=self.user)

            st = self.stat(rfname)
            unique = "%xg%x" % (st.st_dev, st.st_ino)
            modify = time.strftime("%Y%m%d%H%M%S", timefunc(st.st_mtime))

            lines.append(
                f"type={filetype};size={st.st_size};perm={perm};modify={modify};unique={unique}; {instance}")

        line = "\n".join(lines)
        yield line.encode('utf8', self.cmd_channel.unicode_errors)


    def _check_user_status(self):
        "Check if the logged in user is still a valid one. "
        user = self.user['user']
        if user.is_authenticated and user.profile.state in (Profile.BANNED, Profile.SUSPENDED):
            self.cmd_channel.close()
            return


    def extract(self, path_list):
        "Extract relevant info from a path list without breaking on an index error."
        assert isinstance(path_list, list), "path_list should be a list."

        # Project user is under.
        root_project = fetch(len(path_list) >= 1, path_list, 0) or ''

        # Tab picked ( data or results ).
        tab = fetch(len(path_list) >= 2, path_list, 1)

        # Name of specific Data or Job instance.
        name = fetch(len(path_list) >= 3, path_list, 2)

        # Rest of the directory tree to build the actual path with.
        base = None if not len(path_list) >= 3 else path_list[3:]


        return root_project, tab, name, base