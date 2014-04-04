__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError
import os, logging
from optparse import make_option

logger = logging.getLogger("simple-logger")

class Command(BaseCommand):
    help = 'Modifies users'

    option_list = BaseCommand.option_list + (
        make_option('-u', dest='uid',
                    help='Select user by id'),

        make_option('-e', dest='email',
                    help='Select user by email'),

        make_option('-p', dest='passwd',
                    help='sets a new password for the user'),
    )

    def handle(self, *args, **options):
        from biostar.apps.users.models import User
        user = None

        uid = options['uid']
        email = options['email']

        if uid:
             user = User.objects.get(pk=uid)
        elif email:
             user = User.objects.get(email=email)

        passwd = options['passwd']
        if user and passwd:
            set_passwd(user, passwd)

def set_passwd(user, passwd):
    logger.info("settings the password for user %s" % user)
    user.set_password(passwd)
    user.save()

