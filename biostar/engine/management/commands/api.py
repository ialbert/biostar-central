import logging
import hjson
import os
from urllib.parse import urljoin
from urllib.request import urlopen
import requests
import sys

from django.core.management.base import BaseCommand
from django.conf import settings
from django.db.models import Q
from django.shortcuts import reverse
from biostar.engine.models import Analysis

logger = logging.getLogger('engine')

DUMP_DIR = os.path.join(settings.BASE_DIR, "..")
# Override the logger.
logger.setLevel(logging.INFO)


def build_api_url(root_url, uid=None, view="recipe_api_list", api_key=None):

    url = reverse(view, kwargs=dict(uid=uid)) if uid else reverse(view)
    full_url = urljoin(root_url, url) + f"?k={api_key}"

    return full_url


def get_recipes(pid, root_url=None, api_key="", rid=""):
    """
    Return recipes belonging to project 'pid' from api if 'root_url' is given
    else return from database.
    """
    # Filter remote site results by 'pid'
    filter_func = lambda key: recipes[key]["project_uid"] == pid
    # Filter by 'rid' instead if that is given.
    if rid:
        filter_func = lambda key: key == rid

    if root_url:
        # Get the recipes from remote url.
        recipe_api = build_api_url(root_url=root_url, api_key=api_key)
        recipes = hjson.loads(urlopen(url=recipe_api).read().decode())
        # Filter recipes from remote host.
        recipes = list(filter(filter_func, recipes))
    else:
        query = Q(uid=rid, project__uid=pid) if rid else Q(project__uid=pid)
        recipes = Analysis.objects.get_all().filter(query).values_list("uid", flat=True)

    return recipes


def put_recipe(url, files, data, uid=""):

    response = requests.put(url=url, files=files, data=data)
    if response.status_code == 404:
        print(f"*** Recipe id : {uid} does not exist on remote host.")
        sys.exit()
    return response


def recipe_loader(project_dir, api_key="", root_url=None, rid=""):
    """
        Load recipes into api/database from a project found in project_dir.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """
    if not os.path.exists(project_dir):
        print(f"*** Project directory: {project_dir} does not exist.")
        sys.exit()

    # Every subdir in 'project_dir' is a recipe_dir.
    recipe_dirs = [r.name for r in os.scandir(project_dir)]
    # Get the specific recipe to load if given.
    recipe_dirs = list(filter(lambda recipe_uid: recipe_uid == rid, recipe_dirs)) if rid else recipe_dirs

    if not recipe_dirs:
        print(f"*** Project directory :{project_dir} does not contain any recipes.")
        sys.exit()

    def upload(uid, view="recipe_api_template", target_file="", is_json=False, is_image=False):

        target = os.path.join(project_dir, uid, target_file)
        payload = dict(k=api_key)
        mode = "rb" if is_image else "r"
        stream = open(target, mode)
        if root_url:
            # Build api url then send PUT request.
            full_url = build_api_url(root_url=root_url, api_key=api_key, view=view, uid=uid)
            response = put_recipe(url=full_url, files=dict(file=stream), data=payload, uid=uid)
            return response
        # Load data in to the database
        read_file = stream.read()
        if is_image:
            recipe = Analysis.objects.get_all(uid=uid).first()
            open(recipe.image.path, "wb").write(read_file)
        else:
            update_query = dict(json_text=read_file) if is_json else dict(template=read_file)
            Analysis.objects.get_all(uid=uid).update(**update_query)
        return uid

    load_recipe = lambda uid: (upload(uid=uid, target_file="json.hjson", view="recipe_api_json", is_json=True),
                               upload(uid=uid, target_file="image.png", view="recipe_api_image", is_image=True),
                               upload(uid=uid, target_file="template.sh"))
    for recipe_uid in recipe_dirs:
        load_recipe(uid=recipe_uid)
        print(f"Loaded recipe id: {recipe_uid}")

    return recipe_dirs


def recipe_dumper(project_dir, pid, root_url=None, api_key="", rid=""):
    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes to dump into project_dir
    recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key, rid=rid)
    if not recipes:
        msg = f"*** No recipes found for project id={pid}"
        print(msg + (f"and recipe id={rid}." if rid else "."))
        sys.exit()

    def download(uid, is_json=False, view="recipe_api_template", fname="", is_image=False):
        if root_url:
            # Get data from the api url
            fullurl = build_api_url(root_url=root_url, api_key=api_key, view=view, uid=uid)
            data = urlopen(url=fullurl).read()
        else:
            # Get data from database
            recipe = Analysis.objects.get_all(uid=uid).first()
            data = recipe.json_text if is_json else recipe.template
            data = open(recipe.image.path, "rb").read() if is_image else data
        # Make the recipe directory.
        recipe_dir = os.path.join(project_dir, uid)
        os.makedirs(recipe_dir, exist_ok=True)
        # Format data and write to outfile.
        data = hjson.dumps(hjson.loads(data)) if is_json else data
        outfile = os.path.join(recipe_dir, fname)
        mode = "wb" if is_image else "w"
        open(outfile, mode).write(data)
        return outfile

    # Dump json, template, and image for a given recipe
    dump_recipe = lambda uid: (download(uid=uid, fname="json.hjson",  is_json=True, view="recipe_api_json"),
                               download(uid=uid, fname="image.png", is_image=True, view="recipe_api_image"),
                               download(uid=uid, fname="template.sh"))
    for recipe_uid in recipes:
        dump_recipe(uid=recipe_uid)
        print(f"Dumped recipe id: {recipe_uid}")
    return recipes


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
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")
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

        if (root_url and load) and not api_key:
            print("*** --key is required when loading data to remote site.")
            return

        project_dir = os.path.join(root_dir, pid)

        if load:
            loaded = recipe_loader(project_dir=project_dir, root_url=root_url, api_key=api_key, rid=rid)
            view = reverse("project_view", kwargs=dict(uid=pid))
            view = reverse("recipe_view", kwargs=dict(uid=rid)) if rid else view
            print(f"{len(loaded)} recipes loaded into {urljoin(root_url, view) if root_url else 'database'}")
        elif dump:
            dumped = recipe_dumper(root_url=root_url, api_key=api_key, project_dir=project_dir, pid=pid, rid=rid)
            print(f"{len(dumped)} recipes dumped into {project_dir}")

