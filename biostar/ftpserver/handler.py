
import logging
import os
from pyftpdlib.handlers import FTPHandler, DTPHandler
from itertools import chain
from biostar.engine import auth, models


from .util import parse_virtual_path, query_tab


logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)





class BiostarDTPHandler(DTPHandler):
    pass



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

        logger.info(f"file={file}.")
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

        root_project, tab, name, tail = parse_virtual_path(ftppath=path)

        projects = self.fs.projects
        user = self.fs.user["user"]

        # Creating a directory at the root dir
        if root_project and not tab:

            if projects.filter(name=root_project):
                self.respond('550 Directory already exists.')
                return
            else:
                # Create a new project
                project = auth.create_project(user=user, name=root_project)
                self.fs.projects = chain(projects, models.Project.objects.filter(pk=project.pk))

                line = self.fs.fs2ftp(path)
                self.respond('257 "%s" directory created.' % line.replace('"', '""'))

            return path

        if tab == "data" and name and not tail:

            instance = query_tab(tab=tab, project=root_project, name=name, show_instance=True)
            if instance:
                self.respond('550 Directory already exists.')
                return
            else:
                project = self.fs.projects.filter(name=root_project)
                data = auth.create_data(project=project.first(), user=user, name=name)
                self.fs.data = chain(self.fs.data, models.Data.objects.filter(pk=data.pk))

                line = self.fs.fs2ftp(path)
                self.respond('257 "%s" directory created.' % line.replace('"', '""'))

                logger.info(f"path={path}")

            return path

        # Add the data to the tail and update the toc_name.
        logger.info(f"new_dir={path}, path={path}, root_project={root_project}, project={projects.filter(name=root_project)}")

        1/0
        return path


    def ftp_STOR(self, file, mode='w'):

        root_project, tab, name, tail = parse_virtual_path(ftppath=file)

        if tab == 'data' and name and not tail:

            instance = query_tab(tab=tab, project=root_project, name=name, show_instance=True)
            if instance:
                self.respond('550 File already exists.')
                return
            else:

                #data = auth.create_data(project=self.fs.projects.filter(name=root_project).first(),
                #                        name=name)
                #self.fs.data = chain(self.fs.data, models.Data.objects.filter(pk=data.pk))

                testing = super(BiostarFTPHandler, self).ftp_STOR(file=file, mode=mode)
                fd = self.data_channel.file_obj
                logger.info(f"file={file}, testing={testing}, file_obj={fd}")

                return file #os.path.join(data.get_data_dir(), name)

        # Return the real file name here, taken from the name and stuff.
        1/0
