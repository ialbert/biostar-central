import logging
import time
import os

from pyftpdlib.filesystems import AbstractedFS
from django.contrib.auth.models import AnonymousUser
from biostar.engine.models import Data, Job
from biostar.engine import auth
from biostar.accounts.models import Profile

from .util import parse_pk, parse_virtual_path, query_tab, filesystem_mapper


logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)



def fetch_file_info(fname, tab=None, tail=[], pk=None):
    """
    Fetch the actual file path and filetype of a given instance.
    """

    # At /data or /results tab, fname is an ftp path that contains a Data or Job pk.
    pk = pk or parse_pk(fname)

    instance = query_tab(tab=tab, pk=pk, show_instance=True)
    data_file = isinstance(instance, Data) and len(instance.get_files()) <= 1

    suffix = fname if not tail else os.path.join(*tail, fname)
    full_path = os.path.join(instance.get_data_dir(), suffix)
    # There is an extra tail to patch together
    if instance and tail:
        rfname = full_path
        filetype = 'dir' if (os.path.exists(rfname) and os.path.isdir(rfname)) else 'file'

    # No tail to deal with
    elif instance and not tail:
        rfname = instance.get_data_dir()
        is_file = (os.path.isfile(rfname) or data_file or os.path.isfile(full_path))
        filetype = "file" if is_file else "dir"

    else:
        rfname, filetype= fname, "dir"

    return rfname, filetype


def filename(instance, place_holder, tail=[]):

    # Single file data objects handled here
    data_file = isinstance(instance, Data) and len(instance.get_files()) <= 1
    if data_file and instance.get_files()[0]:
        real_file = instance.get_files()[0]
        return real_file if os.path.exists(real_file) else place_holder

    # Directories handled after this point.
    try:
        fname = os.path.join(instance.get_data_dir(), *tail)
        is_file = os.path.exists(fname) and os.path.isfile(fname)
    except Exception as exc:
        logger.error(f"{exc}")
        fname, is_file = None, False

    # Return abspath if its a file, otherwise continue.
    return fname if is_file else place_holder


def dir_list(base_dir, is_instance=False, tail=[]):
    "Return contents of a given base dir ( and a tail)."

    # List the files inside of /data and /results tabs
    if base_dir and is_instance:
        return [os.path.basename(p.path) for p in os.scandir(base_dir)]
    try:
        # List sub-dirs found in /data and /results
        suffix = os.path.join(*tail)
        full_path = os.path.join(base_dir, suffix)
        return [os.path.basename(item.path) for item in os.scandir(full_path)]
    except Exception as exc:
        logger.error(f"{exc}")
        return []


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
        self.timefunc = time.localtime if cmd_channel.use_gmt_times else time.gmtime

        logger.info(f"current_user={current_user}")

        # Dict with all users stored by parent class
        self.user_table = self.cmd_channel.authorizer.user_table

        # Logged in user
        self.username = current_user
        self.user = self.user_table.get(current_user) or dict(user=AnonymousUser)

        # Projects the user has access to
        self.projects = auth.get_project_list(user=self.user["user"])


        # Data and results in project
        self.data = Data.objects.filter(deleted=False, project__in=self.projects,
                                   state__in=(Data.READY, Data.PENDING))

        self.jobs = Job.objects.filter(deleted=False, project__in=self.projects)

        # Get authorizer to set write permission on directories.
        self.authorizer = self.cmd_channel.authorizer


        super(BiostarFileSystem, self).__init__(root, cmd_channel)


    def ftp2fs(self, ftppath):
        """
        :param ftppath: relative virtual path requested by client
        :return:
        """
        abs_ftppath = super(BiostarFileSystem, self).ftp2fs(ftppath)

        logger.info(f"fs={abs_ftppath}, ftppath={ftppath}")

        root_project, tab, pk, tail = parse_virtual_path(ftppath=abs_ftppath)

        if not pk:
            # Only look at actual files, dir are returned as is
            return abs_ftppath

        assert tab in ('data', 'results'), tab
        instance = query_tab(tab=tab, pk=pk, show_instance=True)

        if not instance:
            # We let self.valid_path handle invalid paths
            return abs_ftppath

        # Attempt to return abs file path, otherwise returns place holder.
        return filename(tail=tail, place_holder=abs_ftppath, instance=instance)


    def validpath(self, path):
        """
        Validates a virtual file path when dealing with directory traversal
        and also validates actual file paths when its time to download.

        :param path: current path requested by user
        :return: True if path is valid else False
        """

        logger.info(f"path={path}")

        # Check the user status.
        self._check_user_status()
        real_file, virtual_path = False, False

        # Check if its a virtual path first
        virtual_path = self.validate_virtual_path(path=path)

        # If it's not a virtual ftp path, then it has to be a real file path.
        if not virtual_path:
            real_file =  self.validate_actual_path(path=path)

        return real_file or virtual_path


    def isdir(self, path):


        root_project, tab, pk, tail = parse_virtual_path(ftppath=path)
        is_dir = True

        # If the virtual path does not have a 'tail' then its a dir.
        if not tail or path == self.root:
            return is_dir

        if tail:
            path = os.path.join(query_tab(tab=tab, pk=pk), *tail)
            is_dir = not (os.path.exists(path) and os.path.isfile(path))

        logger.info(f"path={path}, is_dir={is_dir}")
        return is_dir


    def listdir(self, path):
        # This is the root as initialized in the base class

        logger.info(f"path={path}, cwd={self._cwd}, root={self.root}")

        root_project, tab, pk, tail = parse_virtual_path(ftppath=path)

        # List projects when user is at the root
        if path == self.root:
            return filesystem_mapper(queryset=self.projects, tag="Project-")

        # Browse /data or /results inside of a project.
        if self.projects.filter(pk=root_project) and not tab:
            return ['data', 'results']

        # List files in /data or /results
        is_tab = tab and (not pk)
        if self.projects.filter(pk=root_project) and is_tab:
            jobs = self.jobs.filter(project__pk=root_project)
            queryset = self.data.filter(project__pk=root_project) if tab == 'data' else jobs
            return filesystem_mapper(queryset=queryset) or []

        # Take a look at specific instance in /data or /results
        is_instance = pk and (not tail)
        base_dir = query_tab(tab=tab, pk=pk)

        # Dir list returned is different depending on the tab
        return dir_list(base_dir=base_dir, is_instance=is_instance, tail=tail)


    def chdir(self, path):
        """
        Change the current directory.
        """
        # note: process cwd will be reset by the caller

        logger.info(f"path={path}")
        #self._cwd = f"foo({path})
        self._cwd = self.fs2ftp(path)


    def format_mlsx(self, basedir, listing, perms, facts, ignore_err=True):
        logger.info(f"basedir={basedir} listing={listing} facts={facts} perms={perms}")

        lines = []
        root_project, tab, pk, tail = parse_virtual_path(ftppath=basedir)

        is_project = root_project and (not tab)
        perm = self.set_permissions(basedir=basedir)

        for fname in listing:

            if basedir == self.root or is_project:
                pk = root_project if is_project else parse_pk(string=fname)
                project = self.projects.filter(pk=pk).first()
                filetype, rfname = 'dir', project.get_project_dir()
            else:
                rfname, filetype = fetch_file_info(fname=fname, tab=tab, pk=pk, tail=tail)

            st = self.stat(rfname)
            unique = "%xg%x" % (st.st_dev, st.st_ino)
            modify = time.strftime("%Y%m%d%H%M%S", self.timefunc(st.st_mtime))

            lines.append(
                f"type={filetype};size={st.st_size};perm={perm};modify={modify};unique={unique}; {fname}")

        line = "\n".join(lines)
        yield line.encode('utf8', self.cmd_channel.unicode_errors)


    def _check_user_status(self):
        "Check if the logged in user is still valid. "
        user = self.user['user']
        if user.is_authenticated and user.profile.state in (Profile.BANNED, Profile.SUSPENDED):
            self.cmd_channel.close()
            return False


    def validate_virtual_path(self, path):
        "Validate the ftp file path user has requested."

        root_project, tab, pk, tail = parse_virtual_path(ftppath=path)

        # Root project must exist
        if root_project and not self.projects.filter(pk=root_project):
            return False

        # Check user did not leave /data or /results while in project.
        if tab and (tab not in ('data', 'results')):
            return False

        # Data or Job must be valid.
        if pk and not (self.data.filter(pk=pk) or self.jobs.filter(pk=pk)):
            return False

        # The path is valid at this point
        is_valid = True
        if tail:
            # No need to break things while validating the path
            try:
                actual_path = query_tab(tab=tab, pk=pk)
                suffix = os.path.join(*tail)
                is_valid = os.path.exists(os.path.join(actual_path, suffix))
            except Exception as exc:
                logger.error(f"{exc}")
                is_valid = False

        logger.info(f"path={is_valid}")

        return is_valid


    def validate_actual_path(self, path):
        "Make sure user has access to real path requested."


        valid_project_dirs = [p.get_project_dir() for p in self.projects]
        valid_job_dirs = [job.get_data_dir() for job in self.jobs]

        for basedir in valid_project_dirs + valid_job_dirs:
            if path.startswith(basedir):
                logger.info(f"path={path}, basedir={basedir}")
                return True

        return False


    def set_permissions(self, basedir):

        self._check_user_status()

        root_project, tab, pk, tail = parse_virtual_path(ftppath=basedir)


        if tab and (not pk):
            perm = "elrm" if tab == "data" else "elr"

        else:
            perm = "elr"

        if basedir != self.root:
            self.authorizer.override_perm(username=self.username, directory=basedir, perm=perm)


        return perm