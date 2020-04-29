import toml as hjson
import logging
import os

from django.core import management
from django.core.files import File
from django.core.management.base import BaseCommand

from biostar.recipes import auth
from biostar.recipes.models import Project, User
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
        parser.add_argument('--pid',  required=False, help="Project id")

        parser.add_argument('--name', help="Project name")
        parser.add_argument('--demo', action="store_true", default=False, help="Load demo data")

        parser.add_argument('--info', default='', help="File path or text of the project info")

        parser.add_argument('--public', default=False, action="store_true",
                            help="Makes project public")

        parser.add_argument('--update', default=False, action="store_true",
                                        help="Updates the project selected by pid")

    def handle(self, *args, **options):
        pid = options['pid']
        name = options['name'] or "Project Name"
        privacy = Project.PUBLIC if options["public"] else Project.PRIVATE
        update = options["update"]
        info = options["info"]
        demo = options['demo']

        # Set variables needed to execute demo
        if demo:
            update = False
            pid = 'demo'
            privacy = Project.PUBLIC
            info = "Project information goes here."

        # Find project at uid.
        project = Project.objects.filter(uid=pid).first()

        # You can only update existing projects.
        if project and not update:
            error(f"Project with pid={pid} already exists! Set --update to overwrite.")

        # Project description may be a text of a file path.
        if info and os.path.isfile(info):
            # Read the project info from a file.
            info = open(info).read()

        # Select project owner as the first staff user.
        user = User.objects.filter(is_staff=True).first()

        # Create the project.
        pr = auth.create_project(user=user, name=name, uid=pid, text=info, update=update, privacy=privacy)

        if update:
            print(f"*** Updated project uid={pr.uid} name={pr.name}")
        else:
            print(f"*** Created project uid={pr.uid} name={pr.name}")
