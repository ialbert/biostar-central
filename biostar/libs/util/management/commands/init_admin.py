from __future__ import print_function, unicode_literals, absolute_import, division
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
import os

__author__ = 'ialbert'

class Command(BaseCommand):
    help = 'initializes the database with the admin users'

    def handle(self, *args, **options):
        from biostar.apps.accounts.models import User

        #user, flag = User.objects.