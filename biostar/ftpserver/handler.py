
import logging
import time
import os

from pyftpdlib.filesystems import AbstractedFS
from pyftpdlib.handlers import FTPHandler

from django.contrib.auth.models import AnonymousUser
from biostar.engine.models import Project, Data
from biostar.engine import auth

from .authorizer import perm_map

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
        path_list = os.path.split(basedir)


        project = Project.objects.filter(name=path_list[0].replace('/', '')).first()
        data = Data.objects.filter(name= path_list[1].replace('/', ''),
                                   project=project, deleted=False).first()

        suffix = os.path.join(*os.path.split(basedir)[2:], instance)
        full_path = os.path.join(data.get_data_dir(), suffix)
        rfname = list(filter(lambda x : x.endswith(full_path), data.get_files()))[0]

        filetype = 'dir' if rfname and os.path.isdir(rfname) else 'file'


    return project, rfname, filetype



class BiostarFTPHandler(FTPHandler):

    def on_connect(self):
        print("%s:%s connected" % (self.remote_ip, self.remote_port))

    def on_disconnect(self):
        # do something when client disconnects
        pass

    def on_login(self, username):
        # do something when user login

        logger.info(f"user={username}, username={self.username}, auth={self.authenticated}")

        # Tell the filesystem what user is logged in
        # root is the actual directory
        self.fs = self.abstracted_fs(root="/", cmd_channel=self, current_user=username)

    def on_logout(self, username):
        # do something when user logs out

        pass

    def on_file_sent(self, file):
        # do something when a file has been sent
        # Nothing too special here.

        pass

    def on_file_received(self, file):
        # do something when a file has been received


        pass

    def on_incomplete_file_sent(self, file):
        # do something when a file is partially sent
        pass

    def on_incomplete_file_received(self, file):
        # remove partially uploaded files
        pass



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

        path_list = os.path.split(path)
        root_project, current_data = path_list[0], path_list[1]

        project = Project.objects.filter(name=root_project.replace('/', '')).first()

        data = Data.objects.filter(deleted=False, project=project, name=current_data).first()

        if data:
            return [os.path.basename(p) for p in data.get_files()]

        return []



    def chdir(self, path):
        """
        Change the current directory.
        """
        # note: process cwd will be reset by the caller

        logger.info(f"path={path}")
        #self._cwd = f"foo({path})"
        self._cwd = self.fs2ftp(path)


    def fs2ftp(self, fspath):

        logger.info(f"fspath={fspath}")

        return fspath


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


