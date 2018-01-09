import os

from django.core.management.base import BaseCommand

from pyftpdlib.servers import FTPServer
from biostar.engine.handler import EngineFTPHandler, EngineAuthorizer


def get_user_dir(username):

    return



def main():
    # Instantiate a dummy authorizer for managing 'virtual' users
    authorizer = EngineAuthorizer()

    # Define a new user having full r/w permissions and a read-only
    # anonymous user
    authorizer.add_user('user', '12345', './export/media/projects', perm='elradfmwMT')
    authorizer.add_anonymous('./export/media/projects')

    # Instantiate FTP handler class
    handler = EngineFTPHandler
    handler.authorizer = authorizer

    # Define a customized banner (string returned when client connects)
    handler.banner = "Welcome to Biostar-Engine"

    # listen on 0.0.0.0:8080
    address = ('lvh.me', 8080)
    server = FTPServer(address, handler)

    server.max_cons = 256
    server.max_cons_per_ip = 5

    # start ftp server
    server.serve_forever()


class Command(BaseCommand):
    help = 'Manages ftp server'

    def add_arguments(self, parser):
        pass


    def handle(self, *args, **options):
        main()
        pass

