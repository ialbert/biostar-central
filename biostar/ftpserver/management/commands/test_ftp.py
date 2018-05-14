from django.core.management.base import BaseCommand

from ftplib import FTP
from biostar.settings import *


CURRENT_FILE = os.path.abspath(__file__)



class Command(BaseCommand):
    help = "A Test FTP client used to test the server."


    def handle(self, *args, **options):
        import os
        ftp = FTP()

        ftp.connect(host=FTP_HOST, port=FTP_PORT)

        ftp.login(user=ADMINS[0][1], passwd=DEFAULT_ADMIN_PASSWORD)

        # Initial list of projects
        projects = list(ftp.nlst())
        print('\n'.join(projects) + "\n-----")

        project_name = "Test FTP project"
        # Create a new project at the root directory
        ftp.mkd(project_name)

        # Check to see the project has been created
        project_created = project_name in list(ftp.nlst())
        print(f"Project creation and upload successful: {project_created}\n-----")

        # Switch the cwd to the /data directory of newly created project.
        ftp.cwd(os.path.join(project_name, "data"))

        # Creating a dir here will create a Data object in the project.
        ftp.mkd("Test FTP Data")
        ftp.cwd("Test FTP Data")

        # Make another directory
        ftp.mkd("store")
        ftp.cwd("store")

        # Point to the file we want to upload
        file = open(CURRENT_FILE , 'rb')
        fname  = os.path.basename(CURRENT_FILE)
        cmd = f"STOR {fname}"

        # Add directory as a Data object to test project
        ftp.storlines(cmd=cmd, fp=file)

        # Check to see Data object has been created in database
        # and the dir is in the ftp server.
        data_created =  fname in list(ftp.nlst())
        print(f"Data object and files created: {data_created}. Look in {project_name}\n-----")

