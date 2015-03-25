from __future__ import absolute_import, division, print_function, unicode_literals
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from optparse import make_option
import os, logging, glob, json
from biostar3.forum.models import FederatedContent
from biostar3.utils import flatpage

def abspath(*args):
    "Generates absolute paths."
    return os.path.abspath(os.path.join(*args))


class Command(BaseCommand):
    help = '''Adds pages from files in a directory'''

    option_list = BaseCommand.option_list + (
        make_option('--dir', dest="dirname",
                    help='the directory to crawl'),
    )

    def handle(self, *args, **options):
        dirname = options.get("dirname")
        if dirname:
            flatpage.add_all(dirname)

