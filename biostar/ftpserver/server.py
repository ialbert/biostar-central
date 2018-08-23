from pyftpdlib.servers import FTPServer
from pyftpdlib.log import config_logging
from biostar.ftpserver.handler import BiostarFTPHandler
from biostar.ftpserver.authorizer import BiostarAuthorizer
from biostar.ftpserver.filesystem import BiostarFileSystem
from biostar.settings import *
import logging


config_logging(level=logging.DEBUG)

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)


def start():

    # Engine authorization comes from the accounts.
    authorizer = BiostarAuthorizer()
    # Instantiate FTP and DTP handler class.
    handler = BiostarFTPHandler

    handler.authorizer = authorizer
    handler.abstracted_fs = BiostarFileSystem

    # requires SSL for both control and data channel
    #handler.tls_control_required = True
    #handler.tls_data_required = True
    #handler.certfile = 'keycert.pem'

    # Define a customized banner (string returned when client connects)
    handler.banner = "Welcome to Biostar-Engine"

    address = (FTP_HOST, FTP_PORT)
    server = FTPServer(address, handler)

    # FTP connection settings.
    server.max_cons = 256
    server.max_cons_per_ip = 5

    # Start ftp server.
    server.serve_forever()



