import logging
import os
import ftplib

from pyftpdlib.test import unittest
from pyftpdlib.test import MProcessTestFTPd

#from django.test import TestCase
from biostar.accounts.models import User

logger = logging.getLogger("engine")
logger.setLevel(logging.INFO)

USER = 'test@test'
PASSWD = '1234'

#Timeout in seconds
TIMEOUT = 10



class TestFtpFsOperations(unittest.TestCase):

    "test: PWD, CWD, CDUP, SIZE, RNFR, RNTO, DELE, MKD, RMD, MDTM, STAT, MFMT"
    server_class = MProcessTestFTPd
    client_class = ftplib.FTP

    def setUp(self):

        # Create a user
        self.user = User.objects.create_user(username='test', email=USER, password=PASSWD)

        # Start the server
        self.server = self.server_class()
        self.server.start()

        # Start up a client
        self.client = self.client_class(timeout=TIMEOUT)
        self.client.connect(self.server.host, self.server.port)
        self.client.login(USER, PASSWD)

        #
        self.tempfile = os.path.basename(__file__)
        self.tempdir = os.path.basename(tempfile.mkdtemp(dir=HOME))


    def test_mkd(self):
        tempdir = os.path.basename(tempfile.mktemp(dir=HOME))
        dirname = self.client.mkd(tempdir)
        # the 257 response is supposed to include the absolute dirname
        self.assertEqual(dirname, '/' + tempdir)
        # make sure we can't create directories which already exist
        # (probably not really necessary);
        # let's use a try/except statement to avoid leaving behind
        # orphaned temporary directory in the event of a test failure.
        try:
            self.client.mkd(tempdir)
        except ftplib.error_perm:
            os.rmdir(tempdir)  # ok
        else:
            self.fail('ftplib.error_perm not raised.')