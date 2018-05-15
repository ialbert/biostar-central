
import logging
import os
from pyftpdlib.handlers import FTPHandler, TLS_FTPHandler
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

        #TODO: hardcoded "store-"
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

    def make_project_dir(self, root_project, path):
        """
        Create a project in the database and refresh the
        projects tab in the FTP server so it shows up as a directory.
        """
        user = self.fs.user["user"]

        if self.fs.projects.filter(name=root_project):
            self.respond('550 Directory already exists.')
            return
        # Create a new project
        auth.create_project(user=user, name=root_project)

        # Refresh projects tab
        self.fs.projects = auth.get_project_list(user=user)

        line = self.fs.fs2ftp(path)
        self.respond('257 "%s" directory created.' % line.replace('"', '""'))
        logger.info(f"path={path}")
        return path


    def make_data_dir(self, root_project, path, name):

        project = self.fs.projects.filter(name=root_project)

        user = self.fs.user["user"]
        auth.create_data(project=project.first(), user=user, name=name)
        # Refresh /data tab.
        self.fs.data = models.Data.objects.filter(project=project,
                                                  state__in=(models.Data.READY, models.Data.PENDING))
        line = self.fs.fs2ftp(path)
        self.respond('257 "%s" directory created.' % line.replace('"', '""'))
        logger.info(f"path={path}")
        return path


    def ftp_MKD(self, path):

        # Copied from parent class
        # The 257 response is supposed to include the directory
        # name and in case it contains embedded double-quotes
        # they must be doubled (see RFC-959, chapter 7, appendix 2).

        root_project, tab, name, tail = parse_virtual_path(ftppath=path)

        # Creating a project directory at the root dir
        if root_project and not tab:
            return self.make_project_dir(root_project=root_project, path=path)

        # Create a directory inside of the /results or /data tab
        if name:

            instance = query_tab(tab=tab, project=root_project, name=name, show_instance=True)
            if instance and not tail:
                self.respond('550 Directory already exists.')
                return
            if not instance and tab == "results":
                self.respond("550 Only running a recipe will create directories here.")
                return
            elif not instance and tab == "data":
                self.make_data_dir(root_project=root_project, path=path, name=name)
                return

            if instance and tail:
                file = os.path.join(instance.get_data_dir(), *tail)
                if self.is_linked_dir(file=file, data_dir=instance.get_data_dir()):
                    self.respond('550 Can not write to a linked directory.')
                    return

                self.run_as_current_user(self.fs.mkdir, file)
                line = self.fs.fs2ftp(path)
                self.respond('257 "%s" directory created.' % line.replace('"', '""'))

            if tab == "data":
                # Update the toc file and trigger another save.
                instance.make_toc()
                instance.save()
            return path

        else:
            self.respond("550 Can not create a directory here.")
            return


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


    def is_linked_dir(self, file, data_dir):
        "Check if a file in data_dir is inside of a linked directory"

        for fname in os.scandir(data_dir):

            if file.startswith(fname.path) and fname.is_dir():
                if os.path.islink(fname.path):
                    return True
                else:
                    self.is_linked_dir(file=file, data_dir=fname.path)

        return False


    def create_dirs(self, file, instance, tail=[]):
        # Ensure that the sub dirs in tail exist when uploading a file

        not_linked = not self.is_linked_dir(file=file, data_dir=instance.get_data_dir())

        if not_linked:
            file = os.path.join(instance.get_data_dir(), *tail)

            if not os.path.exists(os.path.dirname(file)):
                os.makedirs(os.path.dirname(file), exist_ok=True)

        return file


    def ftp_STOR(self, file, mode='w'):

        root_project, tab, name, tail = parse_virtual_path(ftppath=file)

        if name:

            project = self.fs.projects.filter(name=root_project).first()
            instance = query_tab(tab=tab, project=root_project, name=name, show_instance=True)

            if instance and not tail:
                self.respond('550 File already exists.')
                return

            if instance and tail:
                file = os.path.join(instance.get_data_dir(), *tail)
                # Ensure that the sub dirs in tail exist
                self.create_dirs(file=file, instance=instance, tail=tail)

            elif not instance and tab == "results":
                self.respond('550 Can not upload to the results tab.')
                return

            elif not instance and tab == 'data':
                instance = auth.create_data(project=project, name=name)
                if tail:
                    file = os.path.join(instance.get_data_dir(), *tail)
                    # Ensure that the sub dirs in tail exist
                    self.create_dirs(file=file, instance=instance, tail=tail)
                else:
                    file = os.path.join(instance.get_data_dir(), name)
                # Refresh the data tab
                self.fs.data = models.Data.objects.filter(project=project,
                                                          state__in=(models.Data.READY, models.Data.PENDING))

            if self.is_linked_dir(file=file, data_dir=instance.get_data_dir()):
                self.respond('550 Can not write to a linked directory.')
                return

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



