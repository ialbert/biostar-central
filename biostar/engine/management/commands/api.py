import logging
import hjson
import os
from functools import partial
from urllib.parse import urljoin
from urllib.request import urlopen

from django.core.management.base import BaseCommand

from django.shortcuts import reverse
from biostar.engine.models import Analysis

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)


def get_recipe_api(uid, view="list"):
    view_map = {"list": reverse("recipe_api_list"),
                "json": reverse("recipe_api_json", kwargs=dict(uid=uid)),
                "template": reverse("recipe_api_template", kwargs=dict(uid=uid))}

    return view_map.get("view")


def recipe_loader(root_dir, pid, api_key, root_url=None, rid=""):
    """Load recipes into api/databse from a root_directory."""
    return


def recipe_dumper(target_dir, root_url=None, api_key="", rid=""):
    """
    Dump recipes from the api/database into a target directory.
    """
    def download(uid, outfile="recipe", is_json=False):
        # Get the relative api url.
        url = get_recipe_api(uid=rid, view="json" if is_json else "template")
        # Make the recipe directory.
        recipe_dir = os.path.join(target_dir, uid)
        os.makedirs(recipe_dir, exist_ok=True)
        if root_url:
            fullurl = urljoin(root_url, url) + f"?k={api_key}"
            data = urlopen(url=fullurl).read().decode()
        else:
            recipe = Analysis.objects.get_all(uid=uid)
            data = recipe.json_text if is_json else recipe.template
        # Format data and dump content into file
        data = hjson.dumps(hjson.loads(data)) if is_json else data
        outfile = os.path.join(recipe_dir, outfile)
        open(outfile, "w").write(data)

    download_json = partial(download, outfile="json.hjson", is_json=True)
    download_template = partial(download, outfile="template.sh")
    dump_recipe = lambda uid: (download_template(uid=uid), download_json(uid=uid))
    # Dump into specific recipe dir if set.
    if rid:
        dump_recipe(uid=rid)
        return
    # Scan through the target dir and dump in recipe dirs found within.
    for recipe in os.scandir(target_dir):
        dump_recipe(uid=recipe.name)
    return


class Command(BaseCommand):
    help = 'Dump and load data using api.'

    def add_arguments(self, parser):

        # Load or dump flags
        parser.add_argument('-l', "--load", action="store_true",
                            help="""Load data to url from a directory.
                                    Load recipes to database if --url is not set.""")

        parser.add_argument('-d', "--dump", action="store_true",
                            help="""Dump recipes from a url to directory. 
                                    Dump recipes from database if --url is not set.""")

        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key required to get private projects.")
        parser.add_argument('--dir', default='', help="Directory to store/load data from.")

        parser.add_argument('--pid', type=str, default="all", required=True, help="Project uid to load or dump.")
        parser.add_argument('--rid', type=str, default="", help="Recipe uid to load or dump.")

    def handle(self, *args, **options):

        load = options.get("load")
        dump = options.get("dump")
        root_url = options["url"]
        api_key = options["key"]
        root_dir = options["dir"]
        rid = options["rid"]
        pid = options["pid"]

        if not (load or dump):
            print("*** Set load (-l) or dump (-d) flag.")
            return

        if load and dump:
            print("*** Only one flag can be set.")
            return

        if load:
            recipe_loader(root_dir=root_dir, pid=pid, root_url=None, api_key="", rid="")
        elif dump:
            target_dir = os.path.join(root_dir, pid)
            recipe_dumper(root_url=root_url, api_key=api_key, target_dir=target_dir, rid=rid)


