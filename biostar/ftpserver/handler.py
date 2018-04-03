
import logging
import os
from pyftpdlib.handlers import FTPHandler

from biostar.engine import auth

from .util import parse_virtual_path, filesystem_mapper


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
        # do something after a file has been sent
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


    def ftp_MKD(self, path):

        # Copied from parent class
        # The 257 response is supposed to include the directory
        # name and in case it contains embedded double-quotes
        # they must be doubled (see RFC-959, chapter 7, appendix 2).

        root_project, tab, pk, tail = parse_virtual_path(ftppath=path)

        projects = self.fs.projects
        user = self.fs.user["user"]

        if not root_project:
            # Create a new project
            project = auth.create_project(user=user, name=os.path.basename(path))
            path = filesystem_mapper(instance=project, tag="Project-")

            line = self.fs.fs2ftp(path)
            self.respond('257 "%s" directory created.' % line.replace('"', '""'))

            return self.fs.root

        project = projects.filter(pk=root_project)


        logger.info(f"new_dir={path}, path={path}, root_project={root_project}, projects={projects}")

        1/0

