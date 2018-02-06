from pyftpdlib.servers import FTPServer
from pyftpdlib.log import config_logging
from biostar.ftpserver.handler import BiostarFTPHandler, BiostarAuthorizer, BiostarFileSystem
from biostar.accounts import models
from biostar import  settings
import logging


config_logging(level=logging.DEBUG)

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)


def start():
    # Engine authorization comes from the accounts.
    authorizer = BiostarAuthorizer()

    user = models.User.objects.filter(email="1@lvh.me").first()

    authorizer.add_user("1@lvh.me",settings.DEFAULT_ADMIN_PASSWORD , user=user, perm='elradfmwMT')

    # When parameter user is not specified; we use AnonymousUser
    authorizer.add_user('user', '12345', perm='elradfmwMT')
    authorizer.add_user('user2', '12345', perm='elradfmwMT')

    # Instantiate FTP handler class.
    handler = BiostarFTPHandler
    handler.authorizer = authorizer

    handler.abstracted_fs = BiostarFileSystem

    # Define a customized banner (string returned when client connects)
    handler.banner = "Welcome to Biostar-Engine"

    # Listen on 0.0.0.0:8021
    address = ('lvh.me', 8021)
    server = FTPServer(address, handler)

    # FTP connection settings.
    server.max_cons = 256
    server.max_cons_per_ip = 5

    # Start ftp server.
    server.serve_forever()



