import logging
import time
import os

from pyftpdlib.filesystems import AbstractedFS
from django.contrib.auth.models import AnonymousUser
from biostar.engine.models import Project, Data
from biostar.engine import auth

from .authorizer import perm_map



def split(path):
    path = os.path.normpath(path)
    return [x for x in path.split(os.sep) if x]

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)




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

    if isinstance(instance, str):
        pname, dname = split(basedir)[0], split(basedir)[1]

        suffix = os.path.join(*split(basedir)[2:], instance)
        project = Project.objects.filter(name=pname).first()
        data = Data.objects.filter(name=dname, project=project, deleted=False).first()
        rfname = os.path.join(data.get_data_dir(), suffix)

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

        # Dict with all users
        self.user_table = self.cmd_channel.authorizer.user_table

        # Logged in user
        self.user = self.user_table.get(current_user) or dict(user=AnonymousUser)

        # Projects the user has access to
        self.project_list = auth.get_project_list(user=self.user["user"])

        super(BiostarFileSystem, self).__init__(root, cmd_channel)


    def validpath(self, path):
        logger.info(f"path={path}")

        return True


    def isdir(self, path):
        logger.info(f"path={path}")
        #return True if path in self.project_list else False
        return True


    def listdir(self, path):
        # This is the root as initialized in the base class

        logger.info(f"path={path}, cwd={self._cwd}")
        # List all projects.
        if path == self._root:
            return self.project_list

        # List data belonging to one project
        project = Project.objects.filter(name=path.replace('/', '')).first()
        if project:
            data_list = Data.objects.filter(deleted=False, project=project)
            return data_list

        path_list = split(path)
        root_project, current_data = path_list[0], path_list[1]
        project = Project.objects.filter(name=root_project).first()

        data = Data.objects.filter(deleted=False, project=project, name=current_data).first()

        if data and len(path_list) == 2:
            # Skip the toc when returning dir contents
            return [os.path.basename(p.path) for p in os.scandir(data.get_data_dir())
                    if p.path != data.get_path()]

        suffix = os.path.join(*path_list[2:])
        full_path = os.path.join(data.get_data_dir(), suffix)
        try:
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

        for instance in listing:

            project, rfname, filetype = fetch_file_info(instance, basedir=basedir)

            perm = perm_map(project=project, user=self.user)
            st = self.stat(rfname)
            unique = "%xg%x" % (st.st_dev, st.st_ino)
            modify = time.strftime("%Y%m%d%H%M%S", timefunc(st.st_mtime))

            lines.append(
                f"type={filetype};size={st.st_size};perm={perm};modify={modify};unique={unique}; {instance}")

        line = "\n".join(lines)
        yield line.encode('utf8', self.cmd_channel.unicode_errors)


