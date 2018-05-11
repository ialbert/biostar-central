
import logging
import os
from pyftpdlib.handlers import FTPHandler, DTPHandler
from itertools import chain
from biostar.engine import auth, models


from .util import parse_virtual_path, query_tab


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
        """Only recieve files in data tabs"""

        #TODO: hardcoded "store-1"
        datauid = [s.split("-")[1] for s in file.split(os.sep) if "store-" in s][0]

        data = models.Data.objects.filter(uid=datauid).first()
        if data:
            # Update the toc
            data.make_toc()
            # Trigger another save
            data.save()

        return


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
                auth.create_project(user=user, name=root_project)

                # Refresh projects tab
                self.fs.projects = auth.get_project_list(user=user)

                line = self.fs.fs2ftp(path)
                self.respond('257 "%s" directory created.' % line.replace('"', '""'))

            return path

        if tab == "data" and name:

            instance = query_tab(tab=tab, project=root_project, name=name, show_instance=True)
            project = self.fs.projects.filter(name=root_project)

            if instance and not tail:
                self.respond('550 Directory already exists.')
                return

            if instance and tail:
                file = os.path.join(instance.get_data_dir(), *tail)
                self.run_as_current_user(self.fs.mkdir, file)
                # Remake the toc file
                instance.make_toc()
                instance.save()
                line = self.fs.fs2ftp(file)
                self.respond('257 "%s" directory created.' % line.replace('"', '""'))

            else:
                auth.create_data(project=project.first(), user=user, name=name)
                self.fs.data = models.Data.objects.filter(project=project,
                                                          state__in=(models.Data.READY, models.Data.PENDING))
                line = self.fs.fs2ftp(path)
                self.respond('257 "%s" directory created.' % line.replace('"', '""'))

                logger.info(f"path={path}")

            return path

        # Add the data to the tail and update the toc_name.
        logger.info(f"new_dir={path}, path={path}, root_project={root_project}, project={projects.filter(name=root_project)}")

        1/0
        return path


    def load_dtp(self, file_object, cmd="STOR"):

        if self.data_channel is not None:
            resp = "Data connection already open. Transfer starting."
            self.respond("125 " + resp)
            self.data_channel.file_obj = file_object
            self.data_channel.enable_receiving(self._current_type, cmd)
        else:
            resp = "File status okay. About to open data connection."
            self.respond("150 " + resp)
            self._in_dtp_queue = (file_object, cmd)



    def ftp_STOR(self, file, mode='w'):

        # Make an empty file with the same name in the dest?

        root_project, tab, name, tail = parse_virtual_path(ftppath=file)

        if tab == 'data' and name:

            project = self.fs.projects.filter(name=root_project).first()
            instance = query_tab(tab=tab, project=root_project, name=name, show_instance=True)

            if instance and not tail:
                self.respond('550 File already exists.')
                return

            if instance and tail:
                file = os.path.join(instance.get_data_dir(), *tail)
                # Ensure that the sub dirs in tail exist
                if not os.path.exists(os.path.dirname(file)):
                    os.makedirs(os.path.dirname(file), exist_ok=True)
            else:
                instance = auth.create_data(project=project, name=name)
                file = os.path.join(instance.get_data_dir(), name)
                # Refresh the data tab
                self.fs.data = models.Data.objects.filter(project=project,
                                                          state__in=(models.Data.READY, models.Data.PENDING))
            # Load the stream into the DTP Data Transfer Protocol
            fd = self.run_as_current_user(self.fs.open, file, mode + 'b')
            self.load_dtp(file_object=fd)

            return file

        elif os.path.exists(file):
            self.respond('550 File already exists.')
            return

        else:
            # Can only upload to the /data for now.
            self.respond("550 Can not upload a file here.")
            return



