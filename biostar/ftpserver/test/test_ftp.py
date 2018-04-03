import logging
import os
import ftplib

from pyftpdlib.test import unittest
from pyftpdlib.test import MProcessTestFTPd

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



class TestFtpFsOperations(unittest.TestCase):

    "test: PWD, CWD, CDUP, SIZE, RNFR, RNTO, DELE, MKD, RMD, MDTM, STAT, MFMT"
    server_class = MProcessTestFTPd
    client_class = ftplib.FTP

    def setUp(self):
        self.server = self.server_class()
        self.server.start()
        self.client = self.client_class(timeout=TIMEOUT)
        self.client.connect(self.server.host, self.server.port)
        self.client.login(USER, PASSWD)
        self.tempfile = os.path.basename(touch(TESTFN))
        self.tempdir = os.path.basename(tempfile.mkdtemp(dir=HOME))

    def tearDown(self):
        self.client.close()
        self.server.stop()
        safe_remove(self.tempfile)
        if os.path.exists(self.tempdir):
            shutil.rmtree(self.tempdir)


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