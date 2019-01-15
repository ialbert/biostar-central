
import hjson
from functools import partial
from urllib.parse import urljoin
from urllib.request import urlopen
from django.core.management.base import BaseCommand
from django.conf import settings

from biostar.utils.shortcuts import reverse
import os


def get_base_url():
    return f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}{settings.HTTP_PORT}"


def import_recipes(project_uid, recipe_dict, base_url, base_dir, api_key=""):
    """
    Downloads json and template from API and puts them in project folder.
    """

    def download(url, uid, outfile="recipe", is_json=False):

        # Make the recipe directory
        dir = os.path.join(base_dir, project_uid, uid)
        os.makedirs(dir, exist_ok=True)
        # Get full url and read content
        fullurl = urljoin(base_url, url) + f"?k={api_key}"
        data = urlopen(url=fullurl).read().decode()
        # Format data and dump content into file
        data = hjson.dumps(hjson.loads(data)) if is_json else data
        outfile = os.path.join(dir, outfile)
        open(outfile, "w").write(data)

    for recipe_uid in recipe_dict:
        create_file = partial(download, uid=recipe_uid)
        create_file(outfile="json.hjson", url=recipe_dict[recipe_uid]["json"], is_json=True)
        create_file(outfile="template.sh", url=recipe_dict[recipe_uid]["template"])

    return


class Command(BaseCommand):
    help = 'Import data from api given a base url.'

    def add_arguments(self, parser):

        parser.add_argument('--base', default='', help="Base url to do a reverse look up of api urls.")
        parser.add_argument('--key', default='', help="API key.")
        parser.add_argument('--output', default='', help="Directory to download data to.")

        pass

    def handle(self, *args, **options):

        base_url = options["base"] or get_base_url()
        base_dir = options["output"] or settings.API_DUMP
        api_key = options["key"]
        # Get the project api list
        api_view = reverse("project_api_list")

        # Join the base url with a given api_view and key.
        api_url = urljoin(base_url, api_view) + f"?k={api_key}"

        # Get the json data
        project_json = urlopen(url=api_url).read().decode()
        json_data = hjson.loads(project_json)

        create_files = partial(import_recipes, base_url=base_url, base_dir=base_dir, api_key=api_key)

        # Each key('pid') in json data is a project uid.
        for pid in json_data:
            create_files(project_uid=pid, recipe_dict=json_data[pid]["recipes"])







