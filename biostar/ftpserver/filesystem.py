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


def query_tab(tab, project, name=None, show_instance=False):
    "Return actual path for a path name in /data or /results tab"

    klass_map = {'data': Data, 'results':Job }
    instance = klass_map[tab].objects.filter(deleted=False,
                                             project=project, name=name).first()
    if show_instance:
        return instance

    return None if not instance else instance.get_data_dir()


def fetch_file_info(instance, project=None, tab='data', tail=[], name=None):
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

        suffix = instance if not tail else os.path.join(*tail, instance)
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
        self.projects = auth.get_project_list(user=self.user["user"])

        super(BiostarFileSystem, self).__init__(root, cmd_channel)


    def mkdir(self, path):
        """Create the specified directory."""
        assert isinstance(path, str), path
        logger.info(f"path={path}")

        #os.mkdir(path)


    def ftp2fs(self, ftppath):

        assert isinstance(ftppath, str), ftppath

        if os.path.normpath(self.root) == os.sep:
            fs = os.path.normpath(self.ftpnorm(ftppath))
        else:
            p = self.ftpnorm(ftppath)[1:]
            fs = os.path.normpath(os.path.join(self.root, p))

        logger.info(f"fs={fs}, ftppath={ftppath}")

        path_list = split(fs)
        root_project, tab, name, tail = self.extract(path_list=path_list)

        if not (tail or name):
            # Only look at actual files, dir are returned as is
            return fs

        assert tab in ('data', 'results'), "Can only look at /data or /results of project"

        project = self.projects.filter(uid=root_project).first()
        instance = query_tab(tab=tab, project=project, name=name, show_instance=True)

        if not instance:
            # We let self.valid_path handle invalid paths
            return fs

        # Single file data objects handled here
        data_file = isinstance(instance, Data) and len(instance.get_files()) <= 1
        if data_file and instance.get_files()[0]:
            real_file =  instance.get_files()[0]
            # self.valid_path handles the rest
            return real_file if os.path.exists(real_file) else fs

        # Directories handled after this point.
        base_dir = query_tab(tab=tab, project=project, name=name)

        try:
            fname = os.path.join(base_dir, *tail)
            is_file = os.path.exists(fname) and os.path.isfile(fname)
        except Exception as exc:
            logger.error(f"{exc}")
            fname, is_file = None, False

        # Return abspath if its a file, otherwise continue.
        return fname if is_file else fs


    def validpath(self, path):
        """
        Validates a virtual file path when dealing with directory traversal
        and also validates actual file paths when its time to download.

        :param path: current path requested by user
        :return: True if path is valid else False
        """

        logger.info(f"path={path}")

        # Check the user still has access.
        self._check_user_status()

        path_list = split(path)
        real_file, virtual_path = False, False

        # Check if its a virtual path first
        virtual_path = self.validate_virtual_path(path_list=path_list)

        # If it's not a virtual ftp path, then it has to be a real file path.
        if not virtual_path:
            real_file =  self.validate_actual_path(path=path)

        return real_file or virtual_path


    def isdir(self, path):
        "This actuall matters lol"
        logger.info(f"path={path}")
        return True


    def listdir(self, path):
        # This is the root as initialized in the base class

        logger.info(f"path={path}, cwd={self._cwd}")

        path_list = split(path)
        root_project, tab, name, tail = self.extract(path_list=path_list)

        project = self.projects.filter(uid=root_project).first()
        inside_project = project and len(path_list) == 2

        # List projects when user is at the root
        if path == self.root:
            #TODO: can only space uid out to a max, otherwise it gets truncated.
            return [f"Project-{p.name}{' '*9}({p.uid})" for p in self.projects]

        # We can choose to browse /data or /results inside of a project.
        if project and len(path_list) == 1:
            return ['data', 'results']

        # List project data under /data
        if inside_project and tab == "data":
            #TODO: redo the naming to include uids
            return Data.objects.filter(deleted=False, project=project, 
                                       state__in=(Data.READY, Data.PENDING)) or []

        # List project results under /results
        if inside_project and tab == 'results':
            # TODO: redo the naming to include uids
            return Job.objects.filter(deleted=False, project=project) or []

        # At this point, len(path_list) > 2 and
        # the tab should be in (data, results).
        # The latter is guaranteed by self.valid_path.
        base_dir = query_tab(tab=tab, project=project, name=name)

        # List the files inside of /data and /results tabs
        if base_dir and len(path_list) == 3:
            return [os.path.basename(p.path) for p in os.scandir(base_dir)]

        try:
            # List sub-dirs found in /data and /results
            suffix = os.path.join(*tail)
            full_path = os.path.join(base_dir, suffix)
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

        root_project, tab, name, tail = self.extract(path_list=path_list)

        project = self.projects.filter(uid=root_project).first()

        for instance in listing:
            perm = perm_map(project=instance, user=self.user)
            # Handle parsing project uid from path when path == self.root
            if len(path_list) == 0:
                current = self.projects.filter(uid=self.parse_uid(name=instance))
                
                rfname, filetype = fetch_file_info(current.first(), tab=tab, name=name,
                                                   project=project, tail=tail)
            # Handle /results and /data
            elif len(path_list) == 1:
                filetype, rfname = 'dir', project.get_data_dir()
            else:
                rfname, filetype = fetch_file_info(instance, tab=tab, name=name,
                                                   project=project, tail=tail)

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

    def parse_uid(self, name):
        "Parse a project uid from given project name."
        return name.split(" ")[-1].replace("(","").replace(")","") or None


    def extract(self, path_list):
        "Extract relevant info from a path list without breaking on an index error."
        assert isinstance(path_list, list), "path_list should be a list."

        # Project user is under
        root_project = fetch(len(path_list) >= 1, path_list, 0) or ''

        # Build filesystem using project uid value since names can the same.
        root_project = self.parse_uid(name=root_project)

        # Tab picked ( data or results ).
        tab = fetch(len(path_list) >= 2, path_list, 1)

        # Name of specific Data or Job instance.
        name = fetch(len(path_list) >= 3, path_list, 2)

        # Rest of the directory tree to build the actual path with.
        tail = None if not len(path_list) >= 3 else path_list[3:]

        return root_project, tab, name, tail


    def validate_virtual_path(self, path_list):

        root_project, tab, name, tail = self.extract(path_list=path_list)
        project = self.projects.filter(uid=root_project).first()

        data = Data.objects.filter(name=name, project=project, deleted=False)
        results = Job.objects.filter(name=name, project=project, deleted=False)

        # Root project must exist
        if root_project and not project:
            return False

        # Check user did not leave /data or /results while in project.
        if tab and (tab not in ('data', 'results')):
            return False

        # Data or Job must be valid.
        if name and not (data or results):
            return False

        # The path is valid at this point
        is_valid = True
        if tail:
            # No need to break things while validating the path
            try:
                actual_path = query_tab(tab=tab, project=project, name=name)
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
