import logging
import hjson
import os
import io
from urllib.parse import urljoin
from urllib.request import urlopen
import requests
import sys
from functools import partial

from django.core.management.base import BaseCommand
from django.db.models import Q
from django.shortcuts import reverse
from biostar.engine.models import Analysis, Project, Data
from biostar.engine.api import change_image, get_thumbnail
from biostar.engine import auth
from biostar.accounts.models import User

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)


def build_api_url(root_url, uid=None, view="recipe_api_list", api_key=None):

    url = reverse(view, kwargs=dict(uid=uid)) if uid else reverse(view)
    full_url = urljoin(root_url, url)

    return full_url


def remote_upload(stream, root_url, uid, api_key, view):
    """
    Upload data found in stream to root_url.
    Currently uses PUT requests
    """

    payload = dict(k=api_key)

    # Build api url then send PUT request.
    full_url = build_api_url(root_url=root_url, view=view, uid=uid)
    response = requests.put(url=full_url, files=dict(file=stream), data=payload)
    if response.status_code == 404:
        logger.error(f"*** Object id : {uid} does not exist on remote host.")
        sys.exit()

    return response


def remote_download(root_url, api_key, view, uid, is_image, outfile, is_json):
    """
    Download data found in root_url using GET request.
    """
    mode = "wb" if is_image else "w"
    # Get data from the api url
    fullurl = build_api_url(root_url=root_url, view=view, uid=uid)
    response = requests.get(url=fullurl, params=dict(k=api_key))
    data = response.content if response.status_code == 200 else b""

    # Leave data encoded if its an image
    data = data if is_image else data.decode()
    # Format data and write to outfile.
    data = hjson.dumps(hjson.loads(data)) if is_json else data
    open(outfile, mode).write(data)

    return


def load_db(uid, stream, pid=None, is_json=False, load_recipe=False):
    """
    Load "stream" into database object "uid".
    Loads object as a project by default.
    """
    def project():
        project = Project.objects.get_all(uid=uid).first()
        if not project:
            # Create empty object if not present and populate.
            # Select project owner.
            user = User.objects.filter(is_staff=True).first()
            project = auth.create_project(user=user, name="Project Name", uid=uid)
        conf = hjson.loads(stream.read())
        name = conf.get("settings", {}).get("name", project.name)
        text = conf.get("settings", {}).get("help", project.text)
        Project.objects.get_all(uid=uid).update(name=name, text=text)
        return project

    def recipe():
        recipe = Analysis.objects.get_all(uid=uid).first()
        if not recipe:
            # Create empty object if not present then populate.
            project = Project.objects.get_all(uid=pid).first()
            if not project:
                logger.error(f"*** Project id:{pid} does not exist.")
                logger.error(f"\n*** Run `manage.py api project -load --pid={pid}` to create it.")
                sys.exit()
            recipe = auth.create_analysis(project=project, json_text="", template="", uid=uid, name="Recipe Name")
        if is_json:
            data = hjson.loads(stream.read())
            name = data.get("settings", {}).get("name", recipe.name)
            text = data.get("settings", {}).get("help", recipe.text)
            Analysis.objects.get_all(uid=uid).update(json_text=hjson.dumps(data), name=name, text=text)
        else:
            Analysis.objects.get_all(uid=uid).update(template=stream.read())
        return recipe

    return recipe() if load_recipe else project()


def upload(uid, root_dir, pid=None, root_url=None, api_key="", view="recipe_api_template", fname="",
           is_image=False, load_recipe=False, is_json=False):

    """
    Upload data into a remote host using API.
    Defaults to local database if root_url is None.
    """

    target = os.path.join(root_dir, uid, fname)
    mode = "rb" if is_image else "r"
    if not os.path.exists(target):
        stream = open(get_thumbnail(), mode) if is_image else io.StringIO("")
    else:
        stream = open(target, mode)

    # Upload to remote host when url is set.
    if root_url:
        return remote_upload(stream=stream, root_url=root_url, uid=uid, api_key=api_key, view=view)

    # Update database info

    if is_image:
        # Update image file.
        mtype = Analysis if load_recipe else Project
        obj = mtype.objects.get_all(uid=uid).first()
        return change_image(obj=obj, file_object=stream)

    return load_db(uid=uid, pid=pid, stream=stream, is_json=is_json, load_recipe=load_recipe)


def get_data_placeholder(is_json, is_image, uid):

    if is_image:
        placeholder = open(get_thumbnail(), "rb").read()
    elif is_json:
        data = dict(settings=dict(uid=uid,
                                  name="Object Name",
                                  image=f"{uid}.png",
                                  help="Help Text"))
        placeholder = hjson.dumps(data)

    else:
        placeholder = ""

    return placeholder


def download(uid, root_dir, root_url=None, api_key="", is_json=False, view="recipe_api_template",
             fname="", is_image=False, mtype=Analysis):

    # Get placeholder in case object has no image.
    img_path = lambda o: o.image.path if o.image else get_thumbnail()
    mode = "wb" if is_image else "w"
    # Make output directory.
    outdir = os.path.join(root_dir, uid)
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, fname)

    if root_url:
        remote_download(root_url=root_url, api_key=api_key, view=view, uid=uid,
                        is_image=is_image, outfile=outfile, is_json=is_json)
        return

    # Get data from database
    obj = mtype.objects.get_all(uid=uid).first()
    if not obj:
        data = get_data_placeholder(is_json=is_json, is_image=is_image, uid=uid)
        open(outfile, mode).write(data)
        return

    if is_image:
        data = open(img_path(obj), "rb").read()
    elif is_json:
        data = hjson.dumps(hjson.loads(obj.json_text))
    else:
        data = obj.template

    open(outfile, mode).write(data)
    return outfile


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
        return list(filter(filter_func, recipes))

    query = Q(uid=rid) if rid else Q(project__uid=pid)
    recipes = Analysis.objects.get_all().filter(query)
    if recipes:
        recipes = recipes.values_list("uid", flat=True)
    else:
        # Allows for the creation of 'rid' if it doesn't exist.
        recipes = [rid] if rid else []

    return recipes


def get_image_name(uid, root_url=None, json="conf.hjson", root_dir=None, api_key="", view="recipe_api_json",
                   mtype=Analysis):

    # Get json from url
    if root_url:
        fullurl = build_api_url(root_url=root_url, view=view, uid=uid)
        response = requests.get(url=fullurl, params=dict(k=api_key))
        json_text = response.text if response.status_code == 200 else ""
    # Get json from a file
    elif root_dir:
        path = os.path.join(root_dir, uid, json)
        json_text = open(path).read() if os.path.exists(path) else ""
    # Get json from database
    else:
        obj = mtype.objects.get_all(uid=uid).first()
        json_text = obj.json_text if obj else ""

    json_settings = hjson.loads(json_text).get("settings", {})

    # Get the local image name from "settings" in json.
    # Defaults to uid.png
    name = json_settings.get("image", f"{uid}.png")

    return name


def recipe_loader(root_dir, pid, api_key="", root_url=None, rid=""):
    """
        Load recipes into api/database from a project found in project_dir.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """
    if not os.path.exists(root_dir):
        logger.error(f"*** Project directory: {root_dir} does not exist.")
        sys.exit()

    # Every subdir in 'project_dir' is a recipe_dir.
    recipe_dirs = [r.name for r in os.scandir(root_dir) if r.is_dir()]
    # Get the specific recipe to load if given.
    recipe_dirs = list(filter(lambda recipe_uid: recipe_uid == rid, recipe_dirs)) if rid else recipe_dirs

    # Prepare the main function used to load.
    load = partial(upload, root_dir=root_dir, root_url=root_url, api_key=api_key, pid=pid, load_recipe=True)

    # Get image name from conf file in directory
    img = lambda uid: get_image_name(uid=uid, root_dir=root_dir)
    for recipe_uid in recipe_dirs:
        load(uid=recipe_uid, fname="conf.hjson", view="recipe_api_json", is_json=True)
        load(uid=recipe_uid, fname=img(uid=recipe_uid), view="recipe_api_image", is_image=True)
        load(uid=recipe_uid, fname="template.sh")

        print(f"*** Loaded recipe id: {recipe_uid}")

    return recipe_dirs


def recipe_dumper(root_dir, pid, root_url=None, api_key="", rid=""):
    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes from API or database.
    recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key, rid=rid)

    dump = partial(download, root_url=root_url, root_dir=root_dir, api_key=api_key)

    # Get image name from json on remote host or local database
    img = lambda uid: get_image_name(uid=uid, root_url=root_url, api_key=api_key)
    for recipe_uid in recipes:
        # Dump json, template, and image for a given recipe
        dump(uid=recipe_uid, fname="conf.hjson", is_json=True, view="recipe_api_json")
        dump(uid=recipe_uid, fname=img(uid=recipe_uid), is_image=True, view="recipe_api_image")
        dump(uid=recipe_uid, fname="template.sh")

        print(f"*** Dumped recipe id: {recipe_uid}")
    return recipes


def project_loader(pid, root_dir, root_url=None, api_key=""):
    """
    Load project from root_dir into remote host or local database
    """

    # Prepare function used to upload
    load = partial(upload, uid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key)

    # Get image name from conf file in directory
    img_name = get_image_name(uid=pid, mtype=Project, root_dir=root_dir, json="conf.hjson")

    load(view="project_api_info", fname="conf.hjson")
    load(is_image=True, view="project_api_image", fname=img_name)

    print(f"*** Loaded project ({pid}).")
    return


def project_dumper(pid, root_dir, root_url=None, api_key=""):
    """
    Dump project from remote host or local database into root_dir
    """

    # Prepare function used to download info and images
    dump = partial(download, mtype=Project, uid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key)

    # Get image name from json on remote host or database
    img_name = get_image_name(uid=pid, mtype=Project, root_url=root_url, view="project_api_info")

    # Dump the project json and image
    dump(fname="conf.hjson", view="project_api_info", is_json=True)
    dump(fname=img_name, view="project_api_image", is_image=True)

    print(f"*** Dumped project {pid}: {root_dir}.")
    return


def data_loader(path, pid=None, uid=None, update_toc=False, name="Data Name", type="",
                text=""):
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

    # Slightly different course of action on file and directories.
    isfile = os.path.isfile(path)
    isdir = os.path.isdir(path)

    # The data field is empty.
    if not (isfile or isdir):
        logger.error(f"Path is not a file a directory: {path}")
        return

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
    print(f"*** Creating data: {name}")

    # Create the data.
    auth.create_data(project=project, path=path, type=type, name=name, text=text)

    return


class Command(BaseCommand):
    help = 'Dump and load data using api.'

    def add_api_commands(self, parser):
        """Add default commands to sub commands"""
        parser.add_argument('-l', "--load", action="store_true",
                                                   help="""Load to url from a directory.
                                                        Load to database if --url is not set.""")
        parser.add_argument('-d', "--dump", action="store_true",
                                            help="""Dump from a url to directory. 
                                                    Dump from database if --url is not set.""")
        parser.add_argument("--pid", type=str, default="", help="Project uid to load from or dump to.")
        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")
        parser.add_argument('--dir', default='', help="Directory to store/load data from.")

        return

    def add_arguments(self, parser):
        # Load or dump flags

        subparsers = parser.add_subparsers(help='API sub-commands to choose from.')

        data_parser = subparsers.add_parser("data", help="Load data to local database.")
        data_parser.add_argument("--path", type=str, help="Path to data.")
        data_parser.add_argument("--name", type=str, help="Name of data.")
        data_parser.add_argument("--pid", type=str, default="", help="Project uid to create data in.")
        data_parser.add_argument("--uid", type=str, help="Data uid to load or update.")
        parser.add_argument('--text', default='', help="A file containing the description of the data")
        parser.add_argument('--name', default='', help="Sets the name of the data")
        parser.add_argument('--type', default='data', help="Sets the type of the data")
        data_parser.add_argument("--update_toc", action="store_true", help="Update table of contents for data --uid.")

        recipe_parser = subparsers.add_parser("recipe", help="Load/Dump recipes from/to remote host or local database.")
        recipe_parser.add_argument('--uid', type=str, default="", help="Recipe uid to load or dump.")
        self.add_api_commands(parser=recipe_parser)

        project_parser = subparsers.add_parser("project",
                                               help="Load/Dump projects from/to remote host or local database.")
        self.add_api_commands(parser=project_parser)

    def check_error(self, **options):

        load = options.get("load")
        dump = options.get("dump")
        root_url = options.get("url")
        api_key = options.get("key")
        pid = options.get("pid")
        sys.argv.append("--help")
        msg = None

        if len(sys.argv) == 2 or not pid:
            msg = "[error] --pid needs to be set."

        if not (load or dump):
            msg = "[error] Set load (-l) or dump (-d) flag."

        if load and dump:
            msg = "[error] Only one flag can be set."

        if (root_url and load) and not api_key:
            msg = "[error] --key is required when loading data to remote site."

        if msg:
            logger.error(msg)
            self.run_from_argv(sys.argv)
            sys.exit()

    def handle(self, *args, **options):

        subcommand = sys.argv[2] if len(sys.argv) > 2 else None

        load = options.get("load")
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir") or os.getcwd()
        uid = options.get("uid")
        pid = options.get("pid")
        path = options.get("path")
        name = options.get("name")
        type = options.get("type")
        update_toc = options.get("update_toc")
        project_dir = os.path.join(root_dir, pid)

        if subcommand == "project":
            self.check_error(**options)
            params = dict(pid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key)
            project_loader(**params) if load else project_dumper(**params)

        elif subcommand == "recipe":
            self.check_error(**options)
            params = dict(root_dir=project_dir, root_url=root_url, api_key=api_key, rid=uid, pid=pid)

            recipes = recipe_loader(**params) if load else recipe_dumper(**params)
            self.stdout.write(self.style.SUCCESS(f"{len(recipes)} recipes {'loaded' if load else 'dumped'}."))

        elif subcommand == "data":
            if not (pid or uid):
                logger.error("[error] --pid or --uid need to be set.")
                self.run_from_argv(sys.argv + ["--help"])
                sys.exit()

            data_loader(pid=pid, path=path, uid=uid, update_toc=update_toc, name=name, type=type)
