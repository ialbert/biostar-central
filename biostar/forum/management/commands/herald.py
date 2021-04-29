import logging
from typing import Any
import os, sys
from django.core.management.base import BaseCommand
from biostar.forum.models import Post
from django.conf import settings

logger = logging.getLogger('engine')

LOCK = os.path.join(settings.INDEX_DIR, 'flag')


def publish():
    """
    Publish most recently accepted links
    """
    return


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--publish', action='store_true', default=False,
                            help="Publish most recently accepted herald links")

    def handle(self, *args, **options):
        # Index all un-indexed posts that have a root.
        logger.debug(f"Database: {settings.DATABASE_NAME}")
