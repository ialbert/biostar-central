import os
import requests
from urllib.parse import urljoin
from functools import partial

from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.utils.shortcuts import reverse


def get_base_url():
    return f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}{settings.HTTP_PORT}"


def get_recipe_dirs(base_dir):
    """
    Return list of all recipe dirs given a base dir with structure
    base/project/recipe/
    Returns a Generator Object.
    """
    # Each subdirectory in the project directory is a recipe uid.
    get_dirs = lambda project: [recipe.path for recipe in project]
    # Each subdirectory in the base is a project uid.
    recipe_dirs = (get_dirs(os.scandir(pid.path)) for pid in os.scandir(base_dir))

    # Flatten the nested list of paths.
    return (inner for outter in recipe_dirs for inner in outter)


def export_recipes(recipe_dirs, base_url, api_key=""):
    """
    Export json and template data from list of recipe_dirs to base_url.
    """

    def upload(uid, view, source_name):

        # Build full api url given the view
        full_url = urljoin(base_url, reverse(view, kwargs=dict(uid=uid)))
        # Get the intended file to upload.
        upload_file = dict(file=open(os.path.join(recipe, source_name), "r"))
        # Prepare the payload with the api_key
        payload = dict(k=api_key)
        # Send a PUT request.
        response = requests.put(url=full_url, files=upload_file, data=payload)
        return response

    for recipe in recipe_dirs:
        recipe_uid = os.path.basename(recipe)
        upload(uid=recipe_uid, view="recipe_api_json", source_name="json.hjson")
        upload(uid=recipe_uid, view="recipe_api_template", source_name="template.sh")

    return


class Command(BaseCommand):
    help = 'Export recipe data from base directory to an api using PUT request.'

    def add_arguments(self, parser):

        parser.add_argument('--base', default=get_base_url(),
                            help="Base url to do a reverse look up of api urls. Default: (default: %(default)s)")
        parser.add_argument('--key', default='', help="API key. (default: empty)")
        parser.add_argument('--data', default=settings.API_DUMP, help="""
                                      Base data directory to export project data from.
                                      Its subdirectories are expected to be: /project/recipe.
                                      (default: %(default)s) .""")
        parser.add_argument('--project', default="", help=""" 
                                        Full path to single project directory to crawl and export recipes from.
                                        (default: empty) .""")


        pass

    def handle(self, *args, **options):

        base_url = options["base"]
        base_dir = options["data"]
        api_key = options["key"]
        project_dir = options["project"]

        upload_recipes = partial(export_recipes, base_url=base_url, api_key=api_key)

        # Upload recipes found in a single project and exit out
        if project_dir:
            recipe_dirs = [r.path for r in os.scandir(project_dir) ]
            upload_recipes(recipe_dirs=recipe_dirs)
            return

        # Upload recipes found in multiple projects under the same base_dir
        recipe_dirs = get_recipe_dirs(base_dir=base_dir)
        upload_recipes(recipe_dirs=recipe_dirs)


