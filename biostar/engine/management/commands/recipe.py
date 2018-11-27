import logging
import requests
import io
from urllib.parse import urljoin

import hjson
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.engine import auth
from biostar.utils.shortcuts import reverse
from biostar.engine.models import Project, Analysis

logger = logging.getLogger(settings.LOGGER_NAME)


def update_from_url(root_url, uid, api_key):
    """
    Update the local recipe from remote site.
    """

    params = {"k": api_key}

    # Construct the API urls and join them with the root.

    json_api_url = urljoin(root_url, reverse("api_json", kwargs=dict(uid=uid)))
    template_api_url = urljoin(root_url, reverse("api_template", kwargs=dict(uid=uid)))

    # Get the local recipe
    recipe = Analysis.objects.get_all(uid=uid).first()

    # Get the remote site info
    remote_json = requests.Session().get(json_api_url, data=params)
    remote_template = requests.Session().get(template_api_url, data=params)

    # Update local recipe with remote site info.
    recipe.json_text = remote_json.text
    recipe.template = remote_template.text
    recipe.save()

    #print(local_json, local_template)

    return


def update_to_url(root_url, uid, api_key):
    """
    Update the remote site with local recipe info.
    """

    params = {"k": api_key}

    # Construct the API urls and join them with the root.
    json_api_url = urljoin(root_url, reverse("api_json", kwargs=dict(uid=uid)))
    template_api_url = urljoin(root_url, reverse("api_template", kwargs=dict(uid=uid)))

    # Get the local recipe
    recipe = Analysis.objects.get_all(uid=uid).first()

    # Load the json and template strings as files
    local_json = {'file': io.StringIO(initial_value=recipe.json_text)}
    local_template = {'file': io.StringIO(initial_value=recipe.template)}

    # Use a PUT request to update remote site
    json_request = requests.Session().put(json_api_url, files=local_json, data=params)
    template_request = requests.Session().put(template_api_url, files=local_template, data=params)

    return


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
        parser.add_argument("--image", default='', help="Image path for the recipe.")

        parser.add_argument('--from_url', default='', help="Update local data from remote site.")
        parser.add_argument('--to_url', default='', help="Update remote site with local data.")
        parser.add_argument('--api_key', default='', help="Api key passed to remote site.")

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
        image = options["image"]
        from_url = options["from_url"]
        to_url = options["to_url"]
        api_key = options["api_key"]

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
        # Get the image stream
        stream = open(image, "rb") if image else None

        rec = auth.create_analysis(project=project, json_text=json, template=template, name=name,
                                   text=text, uid=recipe_uid, update=update, stream=stream)

        # Update the local recipe Json and Template from a remote site.
        if from_url:
            update_from_url(root_url=from_url, uid=recipe_uid, api_key=api_key)

        # Update the remote site with local recipe Json and Template.
        elif to_url:
            update_to_url(root_url=to_url, uid=recipe_uid, api_key=api_key)

        if update:
            print(f"*** Updated recipe uid={rec.uid} name={rec.name}")
        else:
            print (f"*** Created recipe uid={rec.uid} name={rec.name}")

