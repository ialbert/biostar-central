from django.core.management.base import BaseCommand

from ftplib import FTP
#from biostar.engine import auth
from biostar.settings import *





class Command(BaseCommand):
    help = "Add users"


    def handle(self, *args, **options):
        import os
        ftp = FTP()

        ftp.connect(host=FTP_HOST, port=FTP_PORT)

        ftp.login(user=ADMINS[0][1], passwd=DEFAULT_ADMIN_PASSWORD)

        #user= User.objects.filter(email=ADMINS[0][1]).first()
        #projects = auth.get_project_list(user=user)

        files = list(ftp.nlst())
        for f in files:
            #print(ftp.cwd(f))
            print(f)

        test_project = files[0]

        # Switch the current working dir to the data directory of test project
        ftp.cwd(os.path.join(test_project, "data"))

        print(ftp.nlst())

        # Try to upload a file.

        file = open(__file__, 'rb')
        fname  = os.path.basename(__file__)
        cmd = f"STOR {fname}"

        ftp.storlines(cmd=cmd, fp=file)

        success = fname in list(ftp.nlst())

        print(f"Upload successful :{success}")






