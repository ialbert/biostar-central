from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
import os
from django.core.mail import send_mail
from biostar.apps.util import mailer
from biostar.server import tasks

class Command(BaseCommand):
    help = 'tests celery tasks'

    def handle(self, *args, **options):
        tasks.test(100, name="Hello!")
        tasks.test.delay(100, name="Hello!")

