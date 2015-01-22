from __future__ import absolute_import, division, print_function, unicode_literals
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from optparse import make_option
import os, logging, glob, json
from biostar3.forum.models import FederatedContent


def abspath(*args):
    "Generates absolute paths."
    return os.path.abspath(os.path.join(*args))


class Command(BaseCommand):
    help = '''Runs quick patches over the data/server. Use it only if you know what you're doing.'''

    DELETE_SQLITE = "delete_sqlite"
    FEDERATED_CONTENT_DIR = "federated_content"

    option_list = BaseCommand.option_list + (
        make_option('--' + DELETE_SQLITE, dest=DELETE_SQLITE, action='store_true', default=False,
                    help='deletes sqlite database'),
        make_option('--' + FEDERATED_CONTENT_DIR, dest=FEDERATED_CONTENT_DIR, default='',
                    help='insert federated content into the database'),
    )

    def handle(self, *args, **options):

        if options[self.DELETE_SQLITE]:
            target = settings.DATABASE_NAME
            print("*** Deleting sqlite database: %s " % target)
            if os.path.isfile(target):
                os.remove(target)
                print("*** Removed: %s" % target)
            else:
                print("*** File not found: %s" % target)

        federated_dir = options[self.FEDERATED_CONTENT_DIR]
        if federated_dir:
            patt = abspath(federated_dir, "*.json")
            for fname in glob.glob(patt):
                content = open(fname).read()

                # Check the the data is valid
                data = json.loads(content)

                domain = data['domain']
                obj_id = data['obj_id']

                print (obj_id)

                obj = FederatedContent(domain=domain, obj_id=obj_id, content=content)
                obj.save()


