from __future__ import absolute_import, division, print_function, unicode_literals
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
import os

class Command(BaseCommand):
    help = 'deletes an sqlite database'

    def handle(self, *args, **options):
        target = settings.DATABASE_NAME
        if os.path.isfile(target):
            os.remove(target)
            print("*** removed: %s" % target)
        else:
            print("*** file not found: %s" % target)