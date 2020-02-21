import logging

from django.core.management.base import BaseCommand
from biostar.forum import markdown

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Used to test markdown rendering'

    def handle(self, *args, **options):


        print(markdown.test())