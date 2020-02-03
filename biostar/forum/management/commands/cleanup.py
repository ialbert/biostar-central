
import logging
from typing import Any

from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from django.conf import settings
from biostar.forum import search

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def handle(self, *args, **options):

        # Index all un-indexed posts that have a root.
        pass