
import logging

from pyftpdlib.authorizers import DummyAuthorizer, AuthenticationFailed
from pyftpdlib.filesystems import AbstractedFS
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.log import config_logging

from django.contrib.auth.models import AnonymousUser
from biostar.accounts import auth as accounts_auth
from biostar.accounts.models import User
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
        self._cwd = root
        self._root = root
        self.cmd_channel = cmd_channel

        logger.info(f"current_user={current_user}")
        # Get current user info
        self.user_table = self.cmd_channel.authorizer.user_table
        self.user = self.user_table.get(current_user) or dict(user=AnonymousUser)
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

        # Return list of public projects when Anonymous user logs in
        return self.project_list

    def chdir(self, path):
        """
        Change the current directory.
        """
        # note: process cwd will be reset by the caller

        logger.info(f"path={path}")
        #self._cwd = f"foo({path})"
        self._cwd = self.fs2ftp(path)

    #def format_list(self, basedir, listing, ignore_err=True):
    #    logger.info(f"listing={listing}")

    def format_mlsx(self, basedir, listing, perms, facts, ignore_err=True):
        logger.info(f"basedir={basedir} listing={listing} facts={facts}")

        lines = []
        for project in listing:

            lines.append(f"type=dir;size=156;perm=r;modify=20071029155301;unique=8012; {project}")

        line = "\n".join(lines)

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


    def validate_authentication(self, username, password, handler):
        """
        Validate user with biostar.accounts and then add to user_table.
        The stored password in the user_table is a hashed-string
        """

        message, valid_django_user = accounts_auth.check_user(email=username, password=password)

        if valid_django_user:

            user = User.objects.filter(email__iexact=username).order_by('-id').first()
            # Store the hash string
            pwd = user.password

            if self.has_user(username=username) and self.user_table[username]['pwd'] == pwd:
                # Check if the user is added to user table with the correct password.
                return

            # Add to user_table since user is a valid biostar.account instance.
            self.add_user(username=username, password=pwd, user=user)
            return

        raise AuthenticationFailed(f"{message}. Resolve issue through website.")


    def get_home_dir(self, username):
        """
        Return the user's home directory.
        Needs to be here because the base class relies on it.
        """
        # Get the project list for the users here and set it as the home_dir
        return '/tmp/this/should/not/be/used/'