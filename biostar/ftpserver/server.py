import os

from django.core.management.base import BaseCommand

from pyftpdlib.servers import FTPServer
from biostar.ftpserver.handler import EngineFTPHandler, EngineAuthorizer, BiostarFileSystem
from biostar.accounts import models
import logging



#logging.basicConfig(level=logging.INFO)

def start():
    # Engine authorization comes from the accounts.
    authorizer = EngineAuthorizer()

    # This is a test user. Will be removed.
    # TODO: iterate and add users from biostar-accounts here
    for user in models.User.objects.all():
        print (user, user.email)

    authorizer.add_user('user', '12345', uid=0, perm='elradfmwMT')

    # Instantiate FTP handler class.
    handler = EngineFTPHandler
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



