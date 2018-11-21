import logging

import hjson
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.engine import auth
from biostar.engine.models import Project, Analysis

logger = logging.getLogger(settings.LOGGER_NAME)


class Bunch():
    def __init__(self, **kwargs):
        self.value = ''
        self.name = self.summary = ''
        self.help = self.type = self.link = ''
        self.__dict__.update(kwargs)


class Command(BaseCommand):
    help = 'Adds recipe to a project (version 2, will replace current)'

    def add_arguments(self, parser):
        parser.add_argument('--id', default='', help="Select project by primary id")
        parser.add_argument('--uid', default='', help="Select project by unique id")
        parser.add_argument('--recipe_uid', required=True, default='', help="Select recipe by unique id")
        parser.add_argument('--json', help="Path to the json to be added")
        parser.add_argument('--template', help="Path to the template to be added")
        parser.add_argument('--text', default='', help="A file containing the description of the data")
        parser.add_argument('--name', default='', help="Sets the name of the data")

        parser.add_argument('--from_url', default='', help="Update local data from remote site.")
        parser.add_argument('--to_url', default='', help="Update remote site with local data.")
        parser.add_argument("--images",  default='', help="Update remote site with local data.")
        parser.add_argument('--update', default=False, action="store_true",
                            help="Updates the project selected by uid")

    def handle(self, *args, **options):

        # Collect the parameters.
        id = options['id']
        uid = options['uid']
        recipe_uid = options["recipe_uid"]
        json = options['json']
        template = options['template']
        text = options['text']
        name = options['name']
        update = options["update"]

        # Select project by id or uid.
        if id:
            project = Project.objects.get_all(id=id).first()
        else:
            project = Project.objects.get_all(uid=uid).first()

        # Get the project.
        recipe = Analysis.objects.get_all(uid=recipe_uid).first()

        # Project must exist if the recipe has not been created yet.
        if project is None and recipe is None:
            logger.error(f"Project not found! id={id} uid={uid}!")
            return

        # Check if recipe exists
        if recipe and not update:
            logger.error(f"Recipe already exists uid={uid}, add the --update flag to update.")
            return

        # Create the recipe
        json = open(json, "r").read() if json else ""
        template = open(template, "r").read() if template else ""
        project = recipe.project if recipe else project

        rec = auth.create_analysis(project=project, json_text=json, template=template, name=name,
                                   text=text, uid=recipe_uid, update=update)

        if update:
            print(f"*** Updated recipe uid={rec.uid} name={rec.name}")
        else:
            print (f"*** Created recipe uid={rec.uid} name={rec.name}")

