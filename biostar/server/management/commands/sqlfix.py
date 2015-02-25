__author__ = 'ialbert'

from django.core.management import call_command
from django.conf import settings
from django.db import connection, transaction
from django.db.models.loading import get_app
from StringIO import StringIO
from django.core.management.base import BaseCommand, CommandError
import os
from optparse import make_option

os.environ['DJANGO_COLORS'] = 'nocolor'

class Command(BaseCommand):
    help = 'Resets the SQL sequence ids'

    option_list = BaseCommand.option_list + (
        make_option('--drop_index', dest='drop_index', action='store_true', default=False, help='drops the index'),
        make_option('--add_index', dest='add_index', action='store_true', default=False, help='adds the index'),
        make_option('--reset', dest='reset', action='store_true', default=False, help='reset sequences'),
    )

    def handle(self, *args, **options):

        if options['drop_index']:
            sql_command("sqldropindexes", "posts")

        if options['add_index']:
            sql_command("sqlindexes", "posts")

        if options['reset']:
            sql_command("sqlsequencereset", "users posts")

def sql_command(command, apps):
    commands = StringIO()
    cursor = connection.cursor()

    targets = apps.split()

    for app in targets:
        label = app.split('.')[-1]
        if get_app(label, emptyOK=True):
            call_command(command, label, stdout=commands)

    sql = commands.getvalue()
    print sql
    cursor.execute(sql)

