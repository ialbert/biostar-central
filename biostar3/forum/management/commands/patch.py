from __future__ import absolute_import, division, print_function, unicode_literals
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from optparse import make_option
import os, logging

class Command(BaseCommand):
    help = '''Runs quick patches over the data/server. Use it only if you know what you're doing.'''

    DELETE_SQLITE = "delete_sqlite"

    option_list = BaseCommand.option_list + (
        make_option('--' + DELETE_SQLITE, dest=DELETE_SQLITE, action='store_true', default=False, help='deletes sqlite database'),

    )

    def handle(self, *args, **options):


        if options[self.DELETE_SQLITE]:
            target = settings.DATABASE_NAME
            print("*** Deleting sqlite database: %s " % target)
            if os.path.isfile(target):
                os.remove(target)
                print("*** Removed: %s" % target)
            else:
                print("*** File not found: %s" % target)