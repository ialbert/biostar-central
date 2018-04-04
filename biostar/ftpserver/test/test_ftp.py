import logging
import os
import ftplib

from pyftpdlib.test import unittest
from pyftpdlib.test import MProcessTestFTPd

from django.test import TestCase

from biostar.ftpserver.util import filesystem_mapper
from biostar.accounts.models import User

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)

USER = 'test@test'
PASSWD = '1234'

#Timeout in seconds
TIMEOUT = 10



class TestFtpFsOperations(TestCase, unittest.TestCase):

    "test: PWD, CWD, CDUP, SIZE, RNFR, RNTO, DELE, MKD, RMD, MDTM, STAT, MFMT"
    server_class = MProcessTestFTPd(addr=('lvh.me', 0000))
    client_class = ftplib.FTP

    def setUp(self):

        # Create a user
        self.user = User.objects.create(username='test', email=USER, password=PASSWD)

        # Start the server
        self.server = self.server_class
        self.server.start()

        # Start up a client
        self.client = self.client_class(timeout=TIMEOUT)
        self.client.connect(self.server.host, self.server.port)
        self.client.login(user=USER, passwd=PASSWD)

        #

    def tearDown(self):
        self.client.close()
        self.server.stop()


    def Xtest_mkd(self):

        temp_project = 'New Project'
        dirname = self.client.mkd(temp_project)
        # the 257 response is supposed to include the absolute dirname
        self.assertTrue(temp_project in dirname)
        1/0
        # make sure we can't create directories which already exist
        # (probably not really necessary);
        # let's use a try/except statement to avoid leaving behind
        # orphaned temporary directory in the event of a test failure.
