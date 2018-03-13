
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


        #self._cwd = root
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

        logger.info(f"path={path}")

        is_project = path  == '/'
        if is_project:
            # Return projects
            return [p.uid for p in self.project_list]

        # Return data contents of a specific project
        project_uid = os.path.split(path)[-1].split('-')[-1]
        project = Project.objects.filter(uid=project_uid).first()
        if not project:
            return  []

        data_list = Data.objects.filter(deleted=False, project=project)

        return [d.uid for d in data_list]


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
        is_project = basedir  == '/'

        if is_project:

            for project_uid in listing:
                project = Project.objects.filter(uid=project_uid).first()
                if not project:
                    continue

                # Virtual fname is what user sees
                # and what we eventually parse in listdir.
                self.add_line(project=project, real_fname=project.get_project_dir(),
                              lines=lines, virtual_fname=f"{project}-{project.uid}")

        else:
            # Show the data as a file for now ( even the directories ).

            for data_uid in listing:
                data = Data.objects.filter(uid=data_uid).first()
                project = data.project
                self.add_line(project=project, real_fname=data.get_data_dir(),
                              lines=lines, virtual_fname=f"{data}-{data.uid}", filetype="file")

        line = "\n".join(lines)
        yield line.encode('utf8', self.cmd_channel.unicode_errors)


    def add_line(self, project, real_fname, virtual_fname, lines=[], filetype="dir"):
        "Append file info to a list as a line"


        if self.cmd_channel.use_gmt_times:
            timefunc = time.gmtime
        else:
            timefunc = time.localtime

        perm = perm_map(project=project, user=self.user)
        st = self.stat(real_fname)

        # Unique fact generated same as pyftpdlib
        unique = "%xg%x" % (st.st_dev, st.st_ino)

        # Show the last time file has been modified.
        modify = time.strftime("%Y%m%d%H%M%S", timefunc(st.st_mtime))

        lines.append(f"type={filetype};size={st.st_size};perm={perm};modify={modify};unique={unique}; {virtual_fname}")
        return


