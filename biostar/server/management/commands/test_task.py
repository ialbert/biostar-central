from __future__ import print_function, unicode_literals, absolute_import, division
import logging

from django.core.management.base import BaseCommand

from biostar.celery import test

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    help = 'tests celery tasks'

    def handle(self, *args, **options):
        logger.info("submitting test task to celery")
        test.delay(100, name="Hello!")


