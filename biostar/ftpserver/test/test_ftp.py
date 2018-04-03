import logging
import os
from ftplib import FTP

from django.test import TestCase


logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)




class CommandsTest(TestCase):


    def setUp(self):
        return


    pass



class DownloadTest(TestCase):

    def setUp(self):
        return




    pass


class UploadTest(TestCase):


    def setUp(self):
        from biostar.accounts.models import User

        self.address = ('lvh.me', 8021)

        self.pswd = "1234"
        self.user = User.objects.create_user(username="test", email="test", password=self.pswd, is_staff=True)


    def test_mkdir(self):

        with FTP(user=self.user.email, passwd=self.pswd, source_address=self.address) as ftp:
            ftp.connect(host=self.address[0], port=self.address[1], timeout=10)
            ftp.login(user=self.user.email, passwd=self.pswd)
            print(ftp.dir())
        1/0
        return





    pass