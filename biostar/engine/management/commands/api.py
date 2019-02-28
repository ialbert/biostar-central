import logging
import hjson
import os
import io
from urllib.parse import urljoin
import requests
import subprocess
import sys

from django.core.management.base import BaseCommand
from django.db.models import Q
from django.shortcuts import reverse
from django.utils import timezone

from django.core.files.base import ContentFile
from biostar.engine.models import Analysis, Project, Data, Job
from biostar.engine.api import change_image, get_thumbnail
from biostar.engine import auth, util
from biostar.accounts.models import User

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)


class Bunch():
    def __init__(self, **kwargs):
        self.value = ''
        self.name = self.summary = ''
        self.help = self.type = self.link = ''
        self.template = self.json_text = ''
        self.text = self.image = ''
        self.json_data = {}
        self.__dict__.update(kwargs)

def build_api_url(root_url, uid=None, view="recipe_api_list", api_key=None):

    url = reverse(view, kwargs=dict(uid=uid)) if uid else reverse(view)
    #TODO Put together params in diffrent way
    full_url = urljoin(root_url, url) + f"?k={api_key}"
    return full_url


def get_recipes(pid, root_url, api_key):
    """
    Return list of recipe uids belonging to project --pid
    """

    if root_url:
        json_url = build_api_url(root_url=root_url, api_key=api_key, view="project_api_info", uid=pid)
        json_text = requests.get(url=json_url, params=dict(k=api_key)).content
        json_data = hjson.loads(json_text)
        # Get the recipes
        recipe_uids = json_data.get("recipes", [])
    else:
        recipes = Analysis.objects.filter(project__uid=pid)
        recipe_uids = recipes.values_list("uid", flat=True)

    if not recipe_uids:
        logger.error(f"No recipes found for pid={pid}")

    return recipe_uids


def generate_fnames(json):

    setting_dict = json.get("settings", {})
    name = setting_dict.get('name')

    name = '_'.join(name.split())
    uid = setting_dict.get("id")
    base = f"{name}-{uid}"

    image = base + '.png'
    template = base + '.sh'
    hjson = base + '.hjson'

    return image, hjson, template


def get_response(root_url, view, uid, api_key=""):
    # Send GET request to view and return a response.
    response = requests.get(url=build_api_url(root_url=root_url, api_key=api_key, uid=uid, view=view))
    if response.status_code != 200:
        logger.error(f"*** Error with remote host with id={uid}: {response.text}")
        sys.exit()
    return response


def put_response(root_url, view, uid, stream, api_key=""):
    # Send PUT request to view with data in stream send as a file.
    full_url = build_api_url(root_url=root_url, view=view, uid=uid, api_key=api_key)
    response = requests.put(url=full_url, files=dict(file=stream), data=dict(k=api_key))
    if response.status_code != 200:
        logger.error(f"*** Error on remote host: {response.text}")
        sys.exit()
    return response


def push_recipe(root_dir,  json_file, api_key="", root_url=None, url_from_json=False):
    """
        Push recipe into api/database from a json file.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """

    # Get conf from json file
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    source = abspath(json_file)
    json_data = hjson.loads(open(source, "r").read())
    image_name, template_name, _ = generate_fnames(json=json_data)

    # Get appropriate parameters from conf
    conf = json_data.get("settings", {})
    rid = conf.get("recipe_uid")
    proj_uid = conf.get("project_uid")
    image = conf.get("image") or image_name
    template = conf.get("template") or template_name
    url = conf.get("url") if url_from_json else root_url
    image_stream = open(abspath(image), "rb")
    template_stream = open(abspath(template), "r")

    if url:
        # Send PUT request to image, json, and template API urls.
        push = lambda view, stream: put_response(root_url=url, view=view, uid=rid, stream=stream, api_key=api_key)
        push(stream=image_stream, view="recipe_api_image")
        push(stream=template_stream, view="recipe_api_template")
        push(stream=open(source, "r"), view="recipe_api_json")
    else:
        project = Project.objects.filter(uid=proj_uid).first()
        recipe = Analysis.objects.filter(uid=rid).first() or Analysis(uid=rid, project=project, owner=project.owner)
        recipe.template = template_stream.read()
        recipe.json_text = hjson.dumps(json_data)
        recipe.name = json_data.get("settings", {}).get("name", recipe.name)
        recipe.text = json_data.get("settings", {}).get("help", recipe.text)
        recipe.image.save(name=image, content=image_stream)
        recipe.save()

    print(f"*** Pushed recipe id={rid} from:{source}")

    return


def push_project(root_dir, json_file, root_url=None, api_key="", url_from_json=False):
    """
    Load projects from root_dir into remote host or local database
    """
    source = os.path.abspath(os.path.join(root_dir, json_file))
    json_data = hjson.loads(open(source, "r").read())
    image_name, _, _ = generate_fnames(json=json_data)

    conf = json_data.get("settings", {})
    image = conf.get('image') or image_name
    image = os.path.abspath(os.path.join(root_dir, image))
    url = conf.get("url") if url_from_json else root_url
    uid = conf.get("uid")
    json_stream = open(source, "r")
    image_stream = open(image, "rb")

    if url:
        # Send PUT request to image, json, and template API urls.
        push = lambda view, stream: put_response(root_url=url, view=view, uid=uid, stream=stream, api_key=api_key)
        push(stream=json_stream, view="project_api_info")
        push(stream=image_stream, view="project_api_image")
    else:
        owner = User.objects.filter(is_superuser=True).first()
        project = Project.objects.filter(uid=uid).first() or Project(uid=uid, owner=owner)
        project.name = conf.get("settings", {}).get("name", project.name)
        project.text = conf.get("settings", {}).get("help", project.text)
        project.image.save(name=image, content=image_stream)
        project.save()

    print(f"*** Pushed project id=({uid}) from:{source}.")


def write_recipe(recipe, root_dir, image):
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))

    # Make output directory.
    os.makedirs(root_dir, exist_ok=True)
    # Write image, template, and json
    fnames = generate_fnames(recipe.json_data)
    img_fname = abspath(fnames[0])
    json_fname = abspath(fnames[1])
    template_fname = abspath(fnames[2])

    open(img_fname, "wb").write(image)
    open(template_fname, "w").write(recipe.template)
    open(json_fname, "w").write(hjson.dumps(recipe.json_data))
    print(f"{json_fname}\n{img_fname}\n{template_fname}")

    return


def write_project(project, root_dir, image):
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))

    os.makedirs(root_dir, exist_ok=True)
    fnames = generate_fnames(project.json_data)
    img_fname = abspath(fnames[0])
    json_fname = abspath(fnames[1])

    # Write image to file
    open(img_fname, "wb").write(image)
    # Write hjson to file.
    open(json_fname, "w").write(hjson.dumps(project.json_data))
    print(f"{json_fname}\n{abspath(img_fname)}")
    return


def pull_recipe(root_dir, rid, url=None, api_key="", save=False):

    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes uid list from API or database.
    get = lambda view: get_response(root_url=url, uid=rid, api_key=api_key, view=view)
    recipe = Analysis.objects.get_all(uid=rid).first()
    image = open(recipe.image.path if recipe.image else get_thumbnail(), "rb").read()

    if url:
        recipe = Bunch(uid=rid)
        json_text = get(view="recipe_api_json").content.decode()
        image = get(view="recipe_api_image").content
        recipe.template = get(view="recipe_api_template").content.decode()
        recipe.json_text = json_text
        recipe.json_data = hjson.loads(json_text)

    if not recipe:
        raise Exception(f"*** Recipe id={rid}does not exist")
    if save:
        print(f"*** Dumping recipe id: {rid}")
        write_recipe(recipe=recipe, root_dir=root_dir, image=image)
    else:
        return recipe

def pull_project(pid, root_dir, url=None, api_key="", save=True):

    """
    Dump project from remote host or local database into root_dir
    """
    project = Project.objects.get_all(uid=pid).first()
    image = open(project.image.path if project.image else get_thumbnail(), "rb").read()

    if url:
        # Get data from remote url.
        project = Bunch(uid=pid)
        image = get_response(root_url=url, api_key=api_key, uid=pid, view="project_api_image").content
        project.json_text = get_response(root_url=url, api_key=api_key, uid=pid, view="project_api_info").content.decode()
        project.json_data = hjson.loads(project.json_text)

    print(f"*** Dumped project {pid}: {root_dir}.")
    if save:
        write_project(project=project, root_dir=root_dir, image=image)
    else:
        return project


def data_loader(path, pid=None, uid=None, update_toc=False, name="Data Name", type="", text=""):
    """
    Load data found in path to database.
    """

    data = Data.objects.get_all(uid=uid).first()
    # Work with existing data.
    if data:
        if update_toc:
            data.make_toc()
            print(f"*** Data id : {uid} table of contents updated.")
        return
    project = Project.objects.get_all(uid=pid).first()
    if not project:
        logger.error(f"Project id: {pid} does not exist.")
        return
    if not path or not os.path.exists(path):
        logger.error(f"--path ({path}) does not exist.")
        return
    # Slightly different course of action on file and directories.
    isdir = os.path.isdir(path)

    # Generate alternate names based on input directory type.
    print(f"*** Project: {project.name} ({project.uid})")
    if isdir:
        print(f"*** Linking directory: {path}")
        altname = os.path.split(path)[0].split(os.sep)[-1]
    else:
        print(f"*** Linking file: {path}")
        altname = os.path.split(path)[-1]

    # Get the text from file
    text = open(text, "r").read() if os.path.exists(text) else ""

    # Select the name.
    name = name or altname

    # Create the data.
    data = auth.create_data(project=project, path=path, type=type, name=name, text=text)
    print(f"*** Created data: name={name}, id={data.uid}")

    return data


def get_json_files(root_dir, json_fname=None):
    """
    Return all existing .hjson or .json files in a directory
    """

    is_json = lambda fname: fname.endswith(".hjson") or fname.endswith(".json")
    json_files = [fname.name for fname in os.scandir(root_dir) if fname.is_file() and is_json(fname.name)]
    if json_fname:
        # Return one json file if provided.
        json_files = [json_fname]

    # Filter for files that exist.
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    recipe_jsons = list(filter(lambda p: os.path.exists(abspath(p)) , json_files))
    return recipe_jsons


def list_obj(mtype="project"):

    mtype_map = dict(project=Project, recipe=Analysis, data=Data, job=Job)

    objs = mtype_map.get(mtype).objects.get_all().order_by("id")[:100]
    for instance in objs:
        print(f'{instance.id}\t{instance.uid}\t{instance.project.uid}\t{instance.name}')

    return


class Command(BaseCommand):
    help = 'Dump and load items using api.'

    def load_data(self, options):
        root_dir = options.get("dir")
        pid = options.get("pid")
        did = options.get("did")
        path = options.get("path") or ""
        name = options.get("name")
        type = options.get("type")
        update_toc = options.get("update_toc")

        if not (pid or did):
            self.stdout.write(self.style.NOTICE("[error] --pid or --did need to be set."))
            self.run_from_argv(sys.argv + ["--help"])
            sys.exit()
        path = path or os.path.abspath(root_dir)
        data = data_loader(pid=pid, path=path, uid=did, update_toc=update_toc, name=name, type=type)
        msg = f"{data.name} loaded into database."
        self.stdout.write(msg=self.style.SUCCESS(msg))
        return

    def manage_push(self, **options):
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir")
        rid = options.get("rid")
        pid = options.get("pid") or ""

        data = options.get("data")
        did = options.get("did")

        json_file = options.get("json")
        url_from_json = options.get("url_from_json")

        # Handle loading data
        if data or did:
            self.load_data(options=options)
            return

        # Get root dir from dir name of json file.
        if json_file:
            full_path = os.path.abspath(json_file)
            root_dir = root_dir or os.path.dirname(full_path)
            json_file = os.path.basename(full_path)

        # Require api key when pushing to remote url
        if (root_url or url_from_json) and not api_key:
            sys.argv.append("--help")
            self.stdout.write(self.style.NOTICE("[error] --key is required when loading data to remote site."))
            self.run_from_argv(sys.argv)
            sys.exit()

        root_dir = root_dir or os.getcwd()
        if not os.path.exists(root_dir):
            logger.error(f"*** Directory: {root_dir} does not exist.")
            sys.exit()

        # Get json files from the root dir.
        json_files = get_json_files(root_dir=root_dir, json_fname=json_file)
        print(f"Pushing files from --dir {root_dir}")
        for fname in json_files:
            json_text = open(os.path.abspath(os.path.join(root_dir, fname)), "r").read()
            recipe_uid = hjson.loads(json_text).get("settings", {}).get("recipe_uid")
            project_uid = hjson.loads(json_text).get("settings", {}).get("project_uid")
            # Skip pushing when rid/pid in json != --rid/--pid given
            if (rid and recipe_uid != rid) or (pid and project_uid != pid):
                continue
            if recipe_uid:
                push_recipe(root_dir=root_dir, root_url=root_url, api_key=api_key, json_file=fname, url_from_json=url_from_json)
                msg = f"Recipe id:{recipe_uid} pushed into {'url' if (url_from_json or root_url) else 'database'}."
            else:
                push_project(root_dir=root_dir, root_url=root_url, api_key=api_key, url_from_json=url_from_json, json_file=fname)
                msg = f"Project id:{recipe_uid} pushed into {'url' if (url_from_json or root_url) else 'database'}."
            logger.info(msg)

        return

    def manage_pull(self, **options):

        pull_recipes = options.get("recipes")
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir") or os.getcwd()
        rid = options.get("rid")
        pid = options.get("pid")

        if (not (pid or rid)) or len(sys.argv) == 3:
            sys.argv.append("--help")
            self.stdout.write(self.style.NOTICE("--pid or --rid is required."))
            self.run_from_argv(sys.argv)
            sys.exit()

        if rid:

            pull_recipe(root_dir=root_dir, url=root_url, api_key=api_key, rid=rid, save=True)

            logger.info(f"Recipe id {rid} dumped into {root_dir}.")
            return

        if pull_recipes:
            # Get recipes belonging to project --pid
            recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key)
            # Get multiple recipes belonging to project --pid
            for uid in recipes:
                pull_recipe(root_dir=root_dir, url=root_url, api_key=api_key, rid=uid, save=True)
                logger.info(f"Recipe id {uid} dumped into {root_dir}.")
            return

        pull_project(pid=pid, root_dir=root_dir, url=root_url, api_key=api_key, save=True)

        logger.info(f"Project id: {pid} dumped into {root_dir}.")
        return

    def add_push_commands(self, parser):

        parser.add_argument('-u', "--url_from_json", action="store_true",
                            help="""Extract url from conf file instead of --url.""")
        parser.add_argument('-d', "--data", action="store_true",
                            help="""Load --dir or --path as a data object of --pid to local database.""")
        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")

        parser.add_argument("--data_from_json", action='store_true', help="Add data found in --json to --pid.")
        parser.add_argument("--jobs", action='store_true', help="Also creates a queued job for the recipe")
        parser.add_argument('--rid', type=str, default="", help="Recipe uid to load.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to load from or dump to.")
        parser.add_argument("--did", type=str, help="Data uid to load or update.")

        parser.add_argument('--dir', default='', help="Base directory to store/load from.")
        parser.add_argument("--path", type=str, help="Path to data.")
        parser.add_argument('--text', default='', help="A file containing the description of the data")
        parser.add_argument('--name', default='', help="Sets the name of the data")
        parser.add_argument('--type', default='data', help="Sets the type of the data")
        parser.add_argument("--update_toc", action="store_true", help="Update table of contents for data --uid.")

        parser.add_argument("--data_root", default="",
                            help="Root directory to data found in conf file when loading project.")
        parser.add_argument('--json', default='', help="""JSON file path relative to --dir to get conf from.""")
        return


    def add_pull_commands(self, parser):
        parser.add_argument('-r', "--recipes", action="store_true",
                            help="""Pull recipes of --pid""")
        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")

        parser.add_argument('--rid', type=str, default="", help="Recipe uid to dump.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to dump.")
        parser.add_argument('--dir', default='', help="Directory to store in.")
        return


    def add_arguments(self, parser):

        subparsers = parser.add_subparsers()

        push_parser = subparsers.add_parser("push", help="""
                                                    Project, Data and Recipe push manager.
                                                    Project:  Create or update project in remote host or database.
                                                    Data:     Create or update data to project --pid in database. 
                                                    Recipe:   Create or update recipe in remote host or database. 
                                                    ."""
                                            )
        self.add_push_commands(parser=push_parser)

        pull_parser = subparsers.add_parser("pull", help="""
                                                    Project and Recipe Job dumper manager.
                                                    Project  : Dump project from remote host or database
                                                    Recipe:  : Dump Recipe from remote host or database.
                                                    """)
        self.add_pull_commands(parser=pull_parser)

    def handle(self, *args, **options):

        subcommand = sys.argv[2] if len(sys.argv) > 2 else None

        listing = options.get("list")

        if len(sys.argv) <= 3:
            sys.argv.append("--help")
            self.run_from_argv(sys.argv)
            sys.exit()

        if listing:
            list_obj(mtype=subcommand)
            return

        if subcommand == "push":
            self.manage_push(**options)
            return

        if subcommand == "pull":
            self.manage_pull(**options)
            return
