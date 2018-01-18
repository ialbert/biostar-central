import logging

from pyftpdlib.authorizers import DummyAuthorizer
from pyftpdlib.filesystems import AbstractedFS
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.log import config_logging

from django.contrib.auth.models import AnonymousUser
from biostar.engine import auth

config_logging(level=logging.DEBUG)

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)



def perm_map():
    return



def project_list(user):

    # Maintain same order as views (doesn't really show)
    projects = auth.get_project_list(user=user).order_by("-sticky", "-privacy")
    projects = projects.order_by("-privacy", "-sticky", "-date", "-id")

    return [p.name for p in projects]


class BiostarFileSystem(AbstractedFS):

    def __init__(self, root, cmd_channel, current_user=None):

        """
         - (str) root: the user "real" home directory (e.g. '/home/user')
         - (instance) cmd_channel: the FTPHandler class instance
         - (str) current_user: username of currently logged in user
        """
        # Set initial current working directory.
        # By default initial cwd is set to "/" to emulate a chroot jail.
        # If a different behavior is desired (e.g. initial cwd = root,
        # to reflect the real filesystem) users overriding this class
        # are responsible to set _cwd attribute as necessary.


        #self._cwd = root
        self._cwd = "/"
        self._root = root
        self.cmd_channel = cmd_channel

        logger.info(f"current_user={current_user}")
        # Get current user info
        self.user_table = self.cmd_channel.authorizer.user_table
        self.user = self.user_table.get(current_user)

        super(BiostarFileSystem, self).__init__(root, cmd_channel)

    def validpath(self, path):
        logger.info(f"path={path}")
        return True

    def isdir(self, path):
        logger.info(f"path={path}")
        return True

    def fs2ftp(self, fspath):
        """Translate a "real" filesystem pathname into equivalent
        absolute "virtual" ftp pathname depending on the user's
        root directory."""
        return fspath


    def ftp2fs(self, ftppath):
        logger.info(f"ftppath={ftppath}")

        #self._cwd = ftppath
        # TODO: the ftppath is going to be

    def listdir(self, path):
        # This is the root as initialized in the base class
        logger.info(f"path={path}")

        if self.user:
            return project_list(user=self.user["user"])

        # Return list of public projects when Anonymous user logs in
        return project_list(user=AnonymousUser)

    def chdir(self, path):
        """
        Change the current directory.
        """
        # note: process cwd will be reset by the caller

        logger.info(f"path={path}")
        #self._cwd = f"foo({path})"
        self._cwd = self.fs2ftp(path)

    def format_list(self, basedir, listing, ignore_err=True):
        logger.info(f"listing={listing}")

    def format_mlsx(self, basedir, listing, perms, facts, ignore_err=True):
        logger.info(f"basedir={basedir} listing={listing} facts={facts}")

        lines = []
        for project in listing:
            #TODO: does it matter if the unique thing is the same for every project
            #TODO: permissons should line up with the access user has to the project
            lines.append(f"type=dir;size=156;perm=r;modify=20071029155301;unique=8012; {project}")

        line = "\n".join(lines)

        #print (line)

        yield line.encode('utf8', self.cmd_channel.unicode_errors)


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
        self.fs = self.abstracted_fs(root="/", cmd_channel=self, current_user=username)

    def on_logout(self, username):
        # do something when user logs out
        pass

    def on_file_sent(self, file):
        # do something when a file has been sent
        #TODO: take basedir as a project and anything in it is a datafile?

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


class BiostarAuthorizer(DummyAuthorizer):
    def add_user(self, username, password, user=AnonymousUser, perm='elr',
                 msg_login="Login successful.", msg_quit="Goodbye."):

        data = {'pwd': str(password),
                'user': user,
                'perm': perm,
                'operms': {},
                'msg_login': str(msg_login),
                'msg_quit': str(msg_quit)
                }
        self.user_table[username] = data

    def get_home_dir(self, username):
        """
        Return the user's home directory.
        Needs to be here because the base class relies on it.
        """
        return '/tmp/this/should/not/be/used/'
