from django.core.management.base import BaseCommand

from ftplib import FTP
#from biostar.engine import auth
from biostar.settings import *





class Command(BaseCommand):
    help = "Add users"


    def handle(self, *args, **options):
        #from biostar.accounts.models import User
        ftp = FTP()     # connect to host, default port

        ftp.connect(host=FTP_HOST, port=FTP_PORT)

        ftp.login(user=ADMINS[0][1], passwd=DEFAULT_ADMIN_PASSWORD)

        #user= User.objects.filter(email=ADMINS[0][1]).first()
        #projects = auth.get_project_list(user=user)

        files = list(ftp.nlst())
        print("ROOT dir")
        for f in files:
            print(ftp.cwd(f))
            print(f)

        #print(files, projects)

