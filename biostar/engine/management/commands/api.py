import logging
import hjson
import os
from functools import partial
from urllib.parse import urljoin
from urllib.request import urlopen

from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.engine.models import Analysis

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)


def get_base_url():
    return f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}{settings.HTTP_PORT}"


def get_api_recipe():

    return


def recipe_loader(**options):

    return


def recipe_dumper(direc, root_url=None, api_key="", recipe=None, project=None):
    """Dump recipes from the api/database into a directory.
    """
    # Get the url

    def download(uid, url=None, outfile="recipe", is_json=False):
        # Make the recipe directory
        recipe_dir = os.path.join(direc, uid)
        os.makedirs(recipe_dir, exist_ok=True)
        if root_url:
            # Get full url and read content
            fullurl = urljoin(root_url, url) + f"?k={api_key}"
            data = urlopen(url=fullurl).read().decode()
        else:
            recipe = Analysis.objects.get_all(uid=uid)
            data = recipe.json_text if is_json else recipe.template
        # Format data and dump content into file
        data = hjson.dumps(hjson.loads(data)) if is_json else data
        outfile = os.path.join(recipe_dir, outfile)
        open(outfile, "w").write(data)

    for recipe_uid in dirs:
        create_file = partial(download, uid=recipe_uid)
        create_file(outfile="json.hjson", url=recipe_dict[recipe_uid]["json"], is_json=True)
        create_file(outfile="template.sh", url=recipe_dict[recipe_uid]["template"])

    return


class Command(BaseCommand):
    help = 'Dump and load data using api.'

    def add_arguments(self, parser):

        # Load or dump flags
        parser.add_argument('-l', "--load", action="store_true",
                            help="Load data to url from a directory.")

        parser.add_argument('-d', "--dump", action="store_true",
                            help="Dump data from a url to directory.")

        parser.add_argument('--url', default="",
                            help="Site url. Dumps to and loads from local database if not provided.")
        parser.add_argument('--key', default='', help="API key.")
        parser.add_argument('--dir', default='', help="Directory to store/load data from.")

        parser.add_argument('--pro', type=str, default="", help="Project uid to load or dump.")
        parser.add_argument('--rec', type=str, default="", help="Recipe uid to load or dump.")

    def handle(self, *args, **options):

        load = options.get("load")
        dump = options.get("dump")
        root_url = options["url"]
        api_key = options["key"]
        direc = options["dir"]
        recipe = options["rec"]
        project = options["pro"]

        if load and dump:
            print("Both load (-l) and dump (-d) flags can not be set.")
            return
        elif not load and not dump:
            print("The load (-l) or dump (-d) flags needs to be set.")
            return
        
        func = recipe_dumper if dump else recipe_loader
        func(root_url=root_url, api_key=api_key, direc=direc, recipe=recipe, project=project)

