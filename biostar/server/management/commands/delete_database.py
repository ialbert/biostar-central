from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
import os


class Command(BaseCommand):
    help = 'deletes an sqlite database'

    def handle(self, *args, **options):
        target = settings.DATABASE_NAME
        if os.path.isfile(target):
            os.remove(target)
        else:
            print("*** file not found: %s" % target)