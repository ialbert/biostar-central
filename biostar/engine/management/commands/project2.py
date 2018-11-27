import hjson
import logging
import os

from django.core import management
from django.core.files import File
from django.core.management.base import BaseCommand

from biostar.engine import auth
from biostar.engine.models import Project, User
import sys


logger = logging.getLogger('engine')


class Bunch():
    def __init__(self, **kwargs):
        self.value = ''
        self.name = self.summary = ''
        self.uid = self.text = ''
        self.user = self.stream = None
        self.__dict__.update(kwargs)


def error(msg):
    logger.error(msg)
    sys.exit()


class Command(BaseCommand):
    help = 'Creates a project.'

    def add_arguments(self, parser):
        parser.add_argument('--uid',  required=True, help="The unique id of the project")

        parser.add_argument('--name', help="The name of the project")

        parser.add_argument('--info', default='', help="File path or text that contains the description")

        parser.add_argument('--public', default=False, action="store_true",
                            help="Makes project public")

        parser.add_argument('--update', default=False, action="store_true",
                        help="Updates the project selected by uid")


    def handle(self, *args, **options):
        uid = options['uid']
        name = options['name'] or "Project Name"
        privacy = Project.PUBLIC if options["public"] else Project.PRIVATE
        update = options["update"]
        info = options["info"]

        # Find project at uid.
        project = Project.objects.get_all(uid=uid).first()

        # You can only update existing projects.
        if project and not update:
            error(f"Project with uid={uid} exists! Set --update to overwrite.")

        # Project description may be a text of a file path.
        if info and os.path.isfile(info):
            # Read the project info from a file.
            info = open(info).read()

        # Select project owner.
        user = User.objects.filter(is_staff=True).first()

        # Create the project.
        pr = auth.create_project(user=user, name=name, uid=uid, text=info, update=update, privacy=privacy)

        if update:
            print(f"*** Updated project uid={pr.uid} name={pr.name}")
        else:
            print (f"*** Created project uid={pr.uid} name={pr.name}")
