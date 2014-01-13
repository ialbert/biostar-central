from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
import os

class Command(BaseCommand):
    help = 'deletes an sqlite database'

    def handle(self, *args, **options):

        os.remove(settings.DATABASE_NAME)