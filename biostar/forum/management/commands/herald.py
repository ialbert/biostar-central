import logging
from typing import Any
import os, sys
from django.template import loader
from django.core.management.base import BaseCommand
from biostar.forum.models import Post, SharedLink
from biostar.forum.util import now
from biostar.accounts.models import User
from biostar.planet.models import BlogPost, Blog
from biostar.forum import auth
from mistune import Markdown
from django.shortcuts import reverse
from django.conf import settings

logger = logging.getLogger('engine')


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--publish', type=int, default=0,
                            help="Create one publication out of N most recently accepted herald submissions")

    def handle(self, *args, **options):
        # Index all un-indexed posts that have a root.
        logger.debug(f"Database: {settings.DATABASE_NAME}")

        publish = options['publish']

        #Herald.objects.update(status=Herald.ACCEPTED)
        #BlogPost.objects.all().delete()

        if publish:
            auth.herald_publisher(limit=publish)
            return
