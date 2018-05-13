from django.core.management.base import BaseCommand

from ftplib import FTP
#from biostar.engine import auth
from biostar.settings import *
from biostar.engine.models import Project

CURRENT_DIR = os.path.join(__file__, '..')


class Command(BaseCommand):
    help = "Add users"


    def handle(self, *args, **options):
        import os
        ftp = FTP()

        ftp.connect(host=FTP_HOST, port=FTP_PORT)

        ftp.login(user=ADMINS[0][1], passwd=DEFAULT_ADMIN_PASSWORD)

        #user= User.objects.filter(email=ADMINS[0][1]).first()
        #projects = auth.get_project_list(user=user)

        # Get initial list of projects
        projects = list(ftp.nlst())
        for f in projects:
            #print(ftp.cwd(f))
            print(f)

        # Test making a project

        test_project = projects[0]

        # Switch the cwd to the /data directory of test_project
        # Creating a dir here will create a Data object in test_project as well.
        ftp.cwd(os.path.join(test_project, "data"))
        print(ftp.nlst())

        # Point to current dir as the target file on local
        file = open(CURRENT_DIR , 'rb')
        fname  = os.path.basename(CURRENT_DIR)
        cmd = f"STOR {fname}"

        # Add file as a Data to the test_project
        ftp.storlines(cmd=cmd, fp=file)

        # Check to see if the ftp server uploaded
        ftp_upload_success = fname in list(ftp.nlst())
        print(f"Upload successful :{ftp_upload_success}")

        # Check to see Data object has been created
        project = Project.objects.filter(name=test_project).first()

        data_obj_created = fname in project.data_set
        print(f"Data object created: {data_obj_created}. Look in {project}")




