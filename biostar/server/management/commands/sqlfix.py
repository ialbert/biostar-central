__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError

class Command(BaseCommand):
    help = 'Resets the SQL sequence ids'

    def handle(self, *args, **options):
        main()


def main():
    commands = StringIO()
    cursor = connection.cursor()
    targets = settings.INSTALLED_APPS

    targets = "biostar.apps.users biostar.apps.posts".split()

    for app in targets:
        label = app.split('.')[-1]
        if get_app(label, emptyOK=True):
            call_command('sqlsequencereset', label, stdout=commands)


    sql = commands.getvalue()
    print sql
    cursor.execute(sql)

