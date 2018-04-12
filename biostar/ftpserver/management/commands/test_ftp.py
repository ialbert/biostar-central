from ftplib import FTP
from biostar.settings import *

ftp = FTP()     # connect to host, default port


ftp.connect(host=FTP_HOST, port=FTP_PORT)


ftp.login()