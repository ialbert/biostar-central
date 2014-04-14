__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError
import os
from optparse import make_option

class Command(BaseCommand):
    help = 'Runs quick patches over the data'

    option_list = BaseCommand.option_list + (
        make_option('--users', dest='users', action='store_true', default=False, help='patches_users'),
    )

    def handle(self, *args, **options):

        if options['users']:
            patch_users()

def patch_users():
    from biostar.apps.users.models import User, Profile
    from biostar.const import DEFAULT_MESSAGES

    users = Profile.objects.all()
    users.update(message_prefs=DEFAULT_MESSAGES)


