import logging
import time
import os

from pyftpdlib.filesystems import AbstractedFS
from django.contrib.auth.models import AnonymousUser
from biostar.engine.models import Project, Data, Job
from biostar.engine import auth
from biostar.accounts.models import Profile

from .authorizer import perm_map



logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)




def index(cond, list, idx):
    # Returns None if cond is not met
    # used to avoid Indexing Errors

    return None if not cond else list[idx]


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


def fetch_file_info(fname, basedir, tab=None, tail=[], pk=None):
    """
    Fetch the actual file path and filetype of a given instance.
    """

    # At /data or /results tab, fname is an ftp path that contains a Data or Job pk.
    if len(split(basedir)) == 2:
        pk = parse_pk(fname)

    instance = query_tab(tab=tab, pk=pk, show_instance=True)
    data_file = isinstance(instance, Data) and len(instance.get_files()) <= 1

    suffix = fname if not tail else os.path.join(*tail, fname)
    # There is an extra tail to patch together
    if instance and tail:
        rfname = os.path.join(instance.get_data_dir(), suffix)
        filetype = 'dir' if (os.path.exists(rfname) and os.path.isdir(rfname)) else 'file'

    # No tail to deal with
    elif instance and not tail:
        rfname = instance.get_data_dir()
        filetype = "file" if data_file else "dir"
    else:
        rfname, filetype= fname, "dir"

    return rfname, filetype


def filesystem_mapper(queryset=None, tag=""):
    "Takes queryset returns correct ftp path to be parsed"

    return [f"{p.pk}.{tag}{p.name}" for p in queryset] or []


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
    root_project = index(len(path_list) >= 1, path_list, 0) or ''

    # Build filesystem using project pk value since names can the same.
    root_project = parse_pk(string=root_project)

    # Tab picked ( data or results ).
    tab = index(len(path_list) >= 2, path_list, 1)

    # Pk of specific Data or Job instance.
    instance = index(len(path_list) >= 3, path_list, 2)

    pk = parse_pk(string=instance)

    # Rest of the directory tree to build the actual path with.
    tail = [] if not len(path_list) >= 3 else path_list[3:]

    return root_project, tab, pk, tail


def get_real_filename(instance, place_holder, tail=[]):

    # Single file data objects handled here
    data_file = isinstance(instance, Data) and len(instance.get_files()) <= 1
    if data_file and instance.get_files()[0]:
        real_file = instance.get_files()[0]
        # self.valid_path handles the rest
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


def get_dir_list(base_dir, is_tab=False, tail=[]):
    "Return contents of a given base dir ( and a tail)."

    # List the files inside of /data and /results tabs
    if base_dir and is_tab:
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

        logger.info(f"current_user={current_user}")

        # Dict with all users stored by parent class
        self.user_table = self.cmd_channel.authorizer.user_table

        # Logged in user
        self.user = self.user_table.get(current_user) or dict(user=AnonymousUser)

        # Projects the user has access to
        self.projects = auth.get_project_list(user=self.user["user"])

        super(BiostarFileSystem, self).__init__(root, cmd_channel)


    def ftp2fs(self, ftppath):

        assert isinstance(ftppath, str), ftppath

        if os.path.normpath(self.root) == os.sep:
            fs = os.path.normpath(self.ftpnorm(ftppath))
        else:
            p = self.ftpnorm(ftppath)[1:]
            fs = os.path.normpath(os.path.join(self.root, p))

        logger.info(f"fs={fs}, ftppath={ftppath}")

        root_project, tab, pk, tail = parse_virtual_path(ftppath=fs)

        if not pk:
            # Only look at actual files, dir are returned as is
            return fs

        assert tab in ('data', 'results'), tab
        instance = query_tab(tab=tab, pk=pk, show_instance=True)

        if not instance:
            # We let self.valid_path handle invalid paths
            return fs

        # Attempt to return abs file path, otherwise returns place holder.
        return get_real_filename(tail=tail, place_holder=fs, instance=instance)


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


        logger.info(f"path={path}")

        root_project, tab, pk, tail = parse_virtual_path(ftppath=path)

        # If the virtual path does not have a 'tail' then its a dir.
        if not tail:
            return True

        return os.path.exists(path) and os.path.isdir(path)


    def listdir(self, path):
        # This is the root as initialized in the base class

        logger.info(f"path={path}, cwd={self._cwd}")

        root_project, tab, pk, tail = parse_virtual_path(ftppath=path)
        path_list = split(path)
        data = Data.objects.filter(deleted=False, project__in=self.projects,
                                   state__in=(Data.READY, Data.PENDING))
        job = Job.objects.filter(deleted=False, project__in=self.projects)

        # List projects when user is at the root
        if path == self.root:
            return filesystem_mapper(queryset=self.projects, tag="Project-")

        # Browse /data or /results inside of a project.
        if self.projects.filter(pk=root_project) and len(path_list) == 1:
            return ['data', 'results']

        # List files in /data or /results
        if self.projects.filter(pk=root_project) and len(path_list) == 2:
            queryset = data if tab == 'data' else job
            return filesystem_mapper(queryset=queryset) or []

        # Take a look at specific instance in /data or /results
        is_tab = len(path_list) == 3
        base_dir = query_tab(tab=tab, pk=pk)

        # Dir list returned is different depending on the tab
        return get_dir_list(base_dir=base_dir, is_tab=is_tab, tail=tail)


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
        timefunc = time.localtime if self.cmd_channel.use_gmt_times else time.gmtime

        root_project, tab, pk, tail = parse_virtual_path(ftppath=basedir)
        perm = perm_map(root_project=root_project, user=self.user)

        for instance in listing:

            if len(split(basedir)) in (0, 1):
                pk = root_project if len(split(basedir)) else parse_pk(string=instance)
                project = self.projects.filter(pk=pk).first()
                filetype, rfname = 'dir', project.get_project_dir()
            else:
                rfname, filetype = fetch_file_info(fname=instance, basedir=basedir, tab=tab, pk=pk, tail=tail)

            st = self.stat(rfname)
            unique = "%xg%x" % (st.st_dev, st.st_ino)
            modify = time.strftime("%Y%m%d%H%M%S", timefunc(st.st_mtime))

            lines.append(
                f"type={filetype};size={st.st_size};perm={perm};modify={modify};unique={unique}; {instance}")

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
        if pk and not (Data.objects.filter(pk=pk, deleted=False) or Job.objects.filter(pk=pk, deleted=False)):
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
        # Jobs found in the projects, that user has access to.
        valid_jobs = Job.objects.filter(deleted=False, project__in=self.projects)

        valid_job_dirs = [job.get_data_dir() for job in valid_jobs]

        for basedir in valid_project_dirs + valid_job_dirs:
            if path.startswith(basedir):
                logger.info(f"path={path}, basedir={basedir}")
                return True

        return False
