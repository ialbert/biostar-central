import logging

from pyftpdlib.authorizers import DummyAuthorizer
from pyftpdlib.filesystems import AbstractedFS
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.log import config_logging

config_logging(level=logging.DEBUG)

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)

def get_user_dir(username):

    return

class BiostarFileSystem(AbstractedFS):

    def __init__(self, root, cmd_channel):
        #TODO: pass user info here ( uid and stuff)

        """
         - (str) root: the user "real" home directory (e.g. '/home/user')
         - (instance) cmd_channel: the FTPHandler class instance
        """
        # Set initial current working directory.
        # By default initial cwd is set to "/" to emulate a chroot jail.
        # If a different behavior is desired (e.g. initial cwd = root,
        # to reflect the real filesystem) users overriding this class
        # are responsible to set _cwd attribute as necessary.
        self._cwd = '/'
        self._root = root
        self.cmd_channel = cmd_channel
        logger.info(f"root={root}")

        super(BiostarFileSystem, self).__init__(root, cmd_channel)

    def validpath(self, path):
        logger.info(f"path={path}")
        return True

    def isdir(self, path):
        logger.info(f"path={path}")
        return True

    def ftp2fs(self, ftppath):
        logger.info(f"ftppath={ftppath}")
        #TODO: take the projects ftp path here and toget path

        #TODO: take out
        self._cwd = ftppath
        print(self._cwd,"SSS")

    def listdir(self, path):
        # This is the root as initialized in the base class
        logger.info(f"path={path}")
        # TODO: list the project here (types will be type=dir)
        # TODO: need to pass the user uid here for that

        return ["music.mp3", "movie.mpg", "image.png"]

    def chdir(self, path):
        """
        Change the current directory.
        """
        # note: process cwd will be reset by the caller
        logger.info(f"path={path}")
        #self._cwd = f"foo({path})"

    def format_list(self, basedir, listing, ignore_err=True):
        logger.info(f"listing={listing}")

    def format_mlsx(self, basedir, listing, perms, facts, ignore_err=True):
        logger.info(f"basedir={basedir} listing={listing}")

        lines = []
        for name in listing:
            lines.append(f"type=file;size=156;perm=r;modify=20071029155301;unique=8012; {name}")

        lines.append(f"type=dir;size=156;perm=r;modify=20071029155301;unique=8012; datadir")
        lines.append(f"type=dir;size=156;perm=r;modify=20071029155301;unique=8012; results")
        lines.append(f"type=dir;size=156;perm=r;modify=20071029155301;unique=8012; nba")
        lines.append(f"type=dir;size=156;perm=r;modify=20071029155301;unique=8012; STUFF")

        line = "\n".join(lines)

        #print (line)

        yield line.encode('utf8', self.cmd_channel.unicode_errors)



class EngineFTPHandler(FTPHandler):
    def on_connect(self):
        print("%s:%s connected" % (self.remote_ip, self.remote_port))

    def on_disconnect(self):
        # do something when client disconnects
        pass

    def on_login(self, username):
        # do something when user login
        pass

    def on_logout(self, username):
        # do something when user logs out
        pass

    def on_file_sent(self, file):
        # do something when a file has been sent
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


class EngineAuthorizer(DummyAuthorizer):
    def add_user(self, username, password, uid, perm='elr',
                 msg_login="Login successful.", msg_quit="Goodbye."):
        data = {'pwd': str(password),
                'uid': uid,
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
        #TODO: should the projects be listed here?
        return '/tmp/this/should/not/be/used/'
