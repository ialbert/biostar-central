import logging
import hjson
import os
from functools import partial
from urllib.parse import urljoin
from urllib.request import urlopen
import requests

from django.core.management.base import BaseCommand

from django.db.models import Q
from django.shortcuts import reverse
from biostar.engine.models import Analysis, Project

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)

DUMP_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "dump")
os.makedirs(DUMP_DIR, exist_ok=True)


def build_api_url(root_url, uid=None, view="recipe_api_list", api_key=None):

    url = reverse(view, kwargs=dict(uid=uid)) if uid else reverse(view)
    full_url = urljoin(root_url, url) + f"?k={api_key}"

    return full_url


def get_recipes(pid, root_url=None, api_key=""):
    """
    Return recipes belonging to project 'pid' from api if 'root_url' is given
    else return from database.
    """
    if root_url:
        # Get the recipes from remote url.
        recipe_api = build_api_url(root_url=root_url, api_key=api_key)
        recipes = hjson.loads(urlopen(url=recipe_api).read().decode())
        # Only get recipes belonging to project "pid".
        recipes = list(filter(lambda key: recipes[key]["project_uid"] == pid, recipes))
    else:
        recipes = Analysis.objects.get_all(project__uid=pid).values_list("uid", flat=True)

    return recipes


def recipe_loader(project_dir, api_key="", root_url=None, rid=""):
    """
        Load recipes into api/database from a project found in project_dir.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """
    # Every subdir in 'project_dir' is a recipe_dir.
    recipe_dirs = [r.name for r in os.scandir(project_dir)]
    # Get a specific recipe to load if given.
    recipe_dirs = filter(lambda recipe_uid: recipe_uid == rid, recipe_dirs) if rid else recipe_dirs

    def upload(uid, view="recipe_api_template", target_file="", is_json=False):

        target = os.path.join(project_dir, uid, target_file)
        # Get the intended file to upload.
        upload_file = dict(file=open(target, "r"))
        payload = dict(k=api_key)
        if root_url:
            # Build api url then send PUT request.
            full_url = build_api_url(root_url=root_url, api_key=api_key, view=view)
            response = requests.put(url=full_url, files=upload_file, data=payload)
            return response

        # Update the recipe.json_text or recipe.template
        update_query = Q(json_text=upload_file["file"].read())
        update_query = update_query if is_json else Q(template=upload_file["file"].read())
        Analysis.objects.get_all(uid=uid).update(update_query)
        return uid

    recipe_loader = lambda uid: (upload(uid=uid, target_file="json.hjson", view="recipe_api_json"),
                                 upload(uid=uid, target_file="template.sh"))
    loaded = list(map(recipe_loader, recipe_dirs))

    return loaded


def recipe_dumper(project_dir, pid, root_url=None, api_key="", rid=""):
    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes to dump into project_dir
    recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key)

    def download(uid, is_json=False, view="recipe_api_template", fname=""):
        # Make the recipe directory.
        recipe_dir = os.path.join(project_dir, uid)
        os.makedirs(recipe_dir, exist_ok=True)
        if root_url:
            # Get data from the api url
            fullurl = build_api_url(root_url=root_url, api_key=api_key, view=view, uid=uid)
            data = urlopen(url=fullurl).read().decode()
        else:
            # Get data from database
            recipe = Analysis.objects.get_all(uid=uid).first()
            data = recipe.json_text if is_json else recipe.template
        # Format data and write to outfile.
        data = hjson.dumps(hjson.loads(data)) if is_json else data
        outfile = os.path.join(recipe_dir, fname)
        open(outfile, "w").write(data)
        return outfile

    # Dump json and template for a given recipe
    dump_recipe = lambda uid: (download(is_json=True, uid=uid, view="recipe_api_json", fname="json.hjson"),
                               download(uid=uid, fname="template.sh"))
    if rid:
        # Dump a specific recipe "rid".
        return dump_recipe(uid=rid)

    dumped = list(map(dump_recipe, recipes))
    return dumped


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

        parser.add_argument('--pid', type=str, required=True, help="Project uid to load or dump.")
        parser.add_argument('--rid', type=str, default="", help="Recipe uid to load or dump.")

    def handle(self, *args, **options):

        load = options.get("load")
        dump = options.get("dump")
        root_url = options["url"]
        api_key = options["key"]
        root_dir = options["dir"] or DUMP_DIR
        rid = options["rid"]
        pid = options["pid"]

        if not (load or dump):
            print("*** Set load (-l) or dump (-d) flag.")
            return

        if load and dump:
            print("*** Only one flag can be set.")
            return

        project_dir = os.path.join(root_dir, pid)
        if load:
            os.makedirs(project_dir, exist_ok=True)
            loaded = recipe_loader(project_dir=project_dir, root_url=root_url, api_key=api_key, rid=rid)
            print(f"{len(loaded)} recipes loaded into {root_url if root_url else 'database'}")
        elif dump:
            dumped = recipe_dumper(root_url=root_url, api_key=api_key, project_dir=project_dir, pid=pid, rid=rid)
            print(f"{len(dumped)} recipes dumped into {project_dir}")

