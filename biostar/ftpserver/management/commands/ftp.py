"""
Takes as input a CVS file with two columns

Email, Name

"""
import csv
import logging
import os

from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.accounts import util
from biostar.accounts.models import User
from biostar.ftpserver import server

logger = logging.getLogger("engine")


class Command(BaseCommand):
    help = "Starts the FTP server"


    def handle(self, *args, **options):
        server.start()

