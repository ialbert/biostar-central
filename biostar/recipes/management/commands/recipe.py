import logging, os
import requests
import io
from urllib.parse import urljoin

import toml as hjson
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.recipes import auth
from biostar.recipes.models import Project, Analysis

logger = logging.getLogger(settings.LOGGER_NAME)


def replace_ext(path, ext='foo'):
    base = os.path.splitext(path)[0]
    newpath = f"{base}.{ext}"
    return newpath

class Command(BaseCommand):
    help = 'Adds recipe to a project (version 2, will replace current)'

    def add_arguments(self, parser):
        parser.add_argument('--pid', default='', required=True, help="Project id.")
        parser.add_argument('--rid', default='', required=True, help="Recipe id.")
        parser.add_argument('--json', help="Recipe json path.")
        parser.add_argument('--template', help="Recipe template path (optional)")
        parser.add_argument('--info', default='', help="Recipe description (optional)")
        parser.add_argument('--name', default='Recipe Name', help="Recipe name")
        parser.add_argument("--image", default='', help="Recipe image path")
        parser.add_argument('--update', default=False, action="store_true",
                            help="Updates the recipe")

    def handle(self, *args, **options):

        # Collect the parameters.
        pid = options['pid']
        rid = options["rid"]
        json = options['json']
        template = options['template'] or replace_ext(json, "sh")
        image = options["image"] or replace_ext(json, "png")
        info = options['info']
        name = options['name']
        update = options["update"]

        # Select project
        project = Project.objects.filter(uid=pid).first() if pid else None

        # Project must exist if the recipe has not been created yet.
        if not project:
            logger.error(f"Project pid={pid} not found!")
            return

        # Get the project.
        recipe = Analysis.objects.filter(uid=rid).first() if rid else None

        # Check if recipe exists
        if recipe and not update:
            logger.error(f"Recipe already exists rid={rid}. Add the --update flag to update.")
            return

        # Create the recipe
        json = open(json, "r").read() if json else ""
        template = open(template, "r").read() if template else ""
        project = recipe.project if recipe else project

        # Get the image stream
        stream = open(image, "rb") if image else None

        rec = auth.create_analysis(project=project, json_text=json, template=template, name=name,
                                   text=info, uid=rid, update=update, stream=stream)

        if update:
            print(f"*** Updated recipe uid={rec.uid} name={rec.name}")
        else:
            print(f"*** Created recipe uid={rec.uid} name={rec.name}")

