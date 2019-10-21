import sys
from urllib.parse import urljoin

import requests
from django.core.management.base import BaseCommand
from django.shortcuts import reverse

from biostar.accounts.models import User
from biostar.recipes import auth
from biostar.recipes.api import *
from biostar.recipes.models import Data

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)


class Bunch():
    def __init__(self, **kwargs):
        self.value = self.uid = ''
        self.name = self.summary = ''
        self.help = self.type = self.link = ''
        self.template = self.json_text = ''
        self.text = self.image = ''
        self.json_data = {}
        self.__dict__.update(kwargs)


def build_api_url(root_url, uid=None, view="recipe_api_list", api_key=None):
    url = reverse(view, kwargs=dict(uid=uid)) if uid else reverse(view)
    # TODO Put together params in diffrent way
    full_url = urljoin(root_url, url) + f"?k={api_key}"
    return full_url


def get_recipes(pid, root_url=None, api_key=""):
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
    name = setting_dict.get('name', 'No name set')

    name = '_'.join(name.split())
    idd = setting_dict.get("id")
    pid = setting_dict.get("project_uid")
    base = f"{name}_{pid}_{idd}"

    image = base + '.png'
    template = base + '.sh'
    hjson = base + '.hjson'

    return image, hjson, template


def get_response(root_url, view, uid, api_key=""):
    # Send GET request to view and return a response.
    url = build_api_url(root_url=root_url, api_key=api_key, uid=uid, view=view)
    response = requests.get(url=url)
    if response.status_code != 200:
        logger.error(f"*** Error with remote host with id={uid}: {response.text}")
        sys.exit()
    return response


def put_response(root_url, view, uid, stream, api_key=""):
    # Send PUT request to view with file in stream sent as a file.
    url = build_api_url(root_url=root_url, view=view, uid=uid, api_key=api_key)
    response = requests.put(url=url, files=dict(file=stream))
    if response.status_code != 200:
        logger.error(f"*** Error on remote host with uid={uid}: {response.text}")
        sys.exit()
    return response


def parse_recipe(json_data):
    """
    Parse recipe information from dict
    """

    image_name,  _, template_name = generate_fnames(json=json_data)

    # Get appropriate parameters
    settings_dict = json_data.get("settings", {})
    rid = settings_dict.get("recipe_uid")
    image = settings_dict.get("image") or image_name
    template = settings_dict.get("template") or template_name
    url = settings_dict.get("url")

    return rid, url, image, template


def parse_project(json_data):
    """
    Parse project information from dict
    """
    image_name, _, _ = generate_fnames(json=json_data)
    pmap = {"private": Project.PRIVATE, "public": Project.PUBLIC}
    settings_dict = json_data.get("settings", {})
    image = settings_dict.get('image') or image_name
    url = settings_dict.get("url")
    pid = settings_dict.get("uid")
    privacy = settings_dict.get("privacy", "").lower() or "private"
    privacy = pmap.get(privacy, Project.PRIVATE)

    return pid, url, privacy, image


def create_recipe(json_file, root_dir, job=False):
    """Create a recipe from a json file."""

    source = os.path.abspath(os.path.join(root_dir, json_file))
    json_data = hjson.loads(open(source, "r").read())
    rid, url, image, template = parse_recipe(json_data)
    recipe = Analysis.objects.filter(uid=rid).first()

    if not recipe:
        # Get an existing project uid
        pid = json_data.get("settings", {}).get("project_uid", )
        name = json_data.get("settings", {}).get("name", "New Recipe")
        text = json_data.get("settings", {}).get("help", "New Text")
        project = Project.objects.filter(uid=pid).first()
        image_stream = open(os.path.abspath(os.path.join(root_dir, image)), "rb")
        template = open(os.path.abspath(os.path.join(root_dir, template)), "r").read()
        recipe = auth.create_analysis(name=name, project=project, text=text, uid=rid,
                                      json_text=hjson.dumps(json_data), template=template, stream=image_stream)
    if job:
        create_job(recipe=recipe)

    return recipe


def create_project(json_file, root_dir):
    """
    Create project from json file.
    """

    source = os.path.abspath(os.path.join(root_dir, json_file))
    json_data = hjson.loads(open(source, "r").read())
    pid, _, privacy, image = parse_project(json_data)

    project = Project.objects.filter(uid=pid).first()

    if not project:
        image_stream = open(os.path.abspath(os.path.join(root_dir, image)), "rb")
        user = User.objects.filter(is_superuser=True).first()
        name = json_data.get("settings", {}).get("name", "New Project")
        text = json_data.get("settings", {}).get("help", "New Text")
        project = auth.create_project(user, name=name, uid=pid, text=text, privacy=privacy)
        project.image.save(name=image, content=image_stream)
        project.save()

    return project


def create_data(root_dir, data_list):
    """Create data from data list, paths are relative to root """

    for data in data_list:
        pid = data.get("project_uid", None)
        # Get data path relative to root_dir
        path = data.get("value", "")
        project = Project.objects.filter(uid=pid).first()
        dtype = data.get("type")
        text = data.get("help")
        name = data.get("name")
        path = path if path.startswith("/") else os.path.join(root_dir, path)

        print(f"*** Creating data in project={pid}, path={path}")
        # Create the data.
        auth.create_data(project=project, path=path, type=dtype, name=name, text=text)

    return


def create_job(recipe):
    """Create a queued job for a recipe"""

    json_data = recipe.json_data

    # Get the data if it exists
    for key, obj in json_data.items():
        if obj.get("source") != "PROJECT":
            continue
        name = obj.get('value', '')
        data = Data.objects.filter(project=recipe.project, name=name).first()
        if not data:
            logger.error(f"Job not created! Missing data:{name} in analysis:{recipe.name}")
            return

        data.fill_dict(obj)

    auth.create_job(analysis=recipe, json_data=json_data)

    return


def push_recipe(root_dir, json_file, api_key="", root_url=None, url_from_json=False):
    """
        Push recipe into api/database from a json file.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """

    # Get conf from json file
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    source = abspath(json_file)
    json_data = hjson.loads(open(source, "r").read())

    rid, url, image, template = parse_recipe(json_data=json_data)
    image_stream = open(abspath(image), "rb")
    template_stream = open(abspath(template), "r")
    url = url if url_from_json else root_url

    if url:
        # Send PUT request to image, json, and template API urls.
        put_response(stream=image_stream, view="recipe_api_image", root_url=url, uid=rid, api_key=api_key)
        put_response(stream=template_stream, view="recipe_api_template", root_url=url, uid=rid, api_key=api_key)
        put_response(stream=open(source, "r"), view="recipe_api_json", root_url=url, uid=rid, api_key=api_key)
    else:
        recipe = Analysis.objects.get(uid=rid)
        recipe.template = template_stream.read()
        recipe.json_text = hjson.dumps(json_data)
        recipe.name = json_data.get("settings", {}).get("name", recipe.name)
        recipe.text = json_data.get("settings", {}).get("help", recipe.text)
        recipe.image.save(name=image, content=image_stream)
        recipe.save()

    logger.info(f"*** Pushed recipe id={rid} from:{source} into {url if url else 'database'}.")

    return


def push_project(root_dir, json_file, root_url=None, api_key="", url_from_json=False):
    """
    Load projects from root_dir into remote host or local database
    """
    source = os.path.abspath(os.path.join(root_dir, json_file))
    json_data = hjson.loads(open(source, "r").read())
    pid, url, _, image = parse_project(json_data)

    url = url if url_from_json else root_url
    json_stream = open(source, "r")
    image_stream = open(os.path.abspath(os.path.join(root_dir, image)), "rb")

    if url:
        # Send PUT request to image, json, and template API urls.
        put_response(stream=json_stream, view="project_api_info", root_url=url, uid=pid, api_key=api_key)
        put_response(stream=image_stream, view="project_api_image", root_url=url, uid=pid, api_key=api_key)
    else:
        project = Project.objects.get(uid=pid)
        project.name = json_data.get("settings", {}).get("name", project.name)
        project.text = json_data.get("settings", {}).get("help", project.text)
        project.image.save(name=image, content=image_stream)
        project.save()

    logger.info(f"*** Pushed project id=({pid}) from:{source} into {url if url else 'database'}.")


def write_recipe(recipe, root_dir, image):

    # Make output directory.
    os.makedirs(root_dir, exist_ok=True)
    # Write image, template, and json
    fnames = generate_fnames(recipe.json_data)
    img_fname = os.path.join(root_dir, fnames[0])
    json_fname = os.path.join(root_dir, fnames[1])
    template_fname = os.path.join(root_dir, fnames[2])

    open(os.path.abspath(img_fname), "wb").write(image)
    open(os.path.abspath(template_fname), "w").write(recipe.template)
    open(os.path.abspath(json_fname), "w").write(hjson.dumps(recipe.json_data))
    print(f"{json_fname}\n{img_fname}\n{template_fname}")

    logger.info(f"Recipe id {recipe.uid} dumped into {root_dir}.")

    return


def write_project(project, root_dir, image):

    os.makedirs(root_dir, exist_ok=True)
    fnames = generate_fnames(project.json_data)
    img_fname = os.path.join(root_dir, fnames[0])
    json_fname = os.path.join(root_dir,fnames[1])

    # Write image to file
    open(os.path.abspath(img_fname), "wb").write(image)
    # Write hjson to file.
    open(os.path.abspath(json_fname), "w").write(hjson.dumps(project.json_data))
    print(f"{json_fname}\n{img_fname}")

    logger.info(f"Project id {project.uid} dumped into {root_dir}.")

    return


def pull_recipe(rid, url=None, api_key=""):
    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes uid list from API or database.
    get = lambda view: get_response(root_url=url, uid=rid, api_key=api_key, view=view)
    recipe = Analysis.objects.filter(uid=rid).first()
    imgpath = recipe.image.path if recipe and recipe.image else get_thumbnail()
    image = open(imgpath, "rb").read()

    if url:
        recipe = Bunch(uid=rid)
        json_text = get(view="recipe_api_json").content.decode()
        image = get(view="recipe_api_image").content
        recipe.template = get(view="recipe_api_template").content.decode()
        recipe.json_text = json_text
        recipe.json_data = hjson.loads(json_text)

    return recipe, image


def pull_project(pid, url=None, api_key=""):
    """
    Dump project from remote host or local database into root_dir
    """
    project = Project.objects.filter(uid=pid).first()
    imgpath = project.image.path if project and project.image else get_thumbnail()
    image = open(imgpath, "rb").read()

    if url:
        # Get data from remote url.
        project = Bunch(uid=pid)
        image = get_response(root_url=url, api_key=api_key, uid=pid, view="project_api_image").content
        project.json_text = get_response(root_url=url, api_key=api_key, uid=pid,
                                         view="project_api_info").content.decode()
        project.json_data = hjson.loads(project.json_text)

    return project, image


def get_json_files(root_dir, json_fname=None):
    """
    Return all existing .hjson or .json files in a directory
    """

    is_json = lambda fname: fname.is_file() and fname.name.endswith(".hjson") or fname.name.endswith(".json")
    json_files = [fname.name for fname in os.scandir(root_dir) if is_json(fname)]
    if json_fname:
        # Return one json file if provided.
        json_files = [json_fname]

    # Filter for files that exist.
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    recipe_jsons = list(filter(lambda p: os.path.exists(abspath(p)), json_files))
    return recipe_jsons


def list_ids(url=None, api_key=""):
    """
    Returns a listing of all project and recipe ids.
    """
    if url:
        url = build_api_url(root_url=url, api_key=api_key, view="api_list")
        response = requests.get(url=url)
        text = response.text
    else:
        text = tabular_list()

    print(text)


class Command(BaseCommand):
    help = 'Dump and load items using api.'

    def manage_push(self, **options):
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir")
        rid = options.get("rid")
        pid = options.get("pid") or ""

        json_file = options.get("json")
        url_from_json = options.get("url_from_json")

        # Get root dir from directory name with json file.
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
        for fname in json_files:
            json_text = open(os.path.abspath(os.path.join(root_dir, fname)), "r").read()
            recipe_uid = hjson.loads(json_text).get("settings", {}).get("recipe_uid")
            project_uid = hjson.loads(json_text).get("settings", {}).get("project_uid")
            # Skip pushing when rid/pid in json != --rid/--pid given
            if (rid and recipe_uid != rid) or (pid and project_uid != pid):
                continue
            if recipe_uid:
                push_recipe(root_dir=root_dir, root_url=root_url, api_key=api_key, json_file=fname,
                            url_from_json=url_from_json)
            else:
                push_project(root_dir=root_dir, root_url=root_url, api_key=api_key, url_from_json=url_from_json,
                             json_file=fname)
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
            recipe, image = pull_recipe(url=root_url, api_key=api_key, rid=rid)
            if recipe:
                write_recipe(recipe=recipe, root_dir=root_dir, image=image)
            return

        if pull_recipes:
            # Get recipes belonging to project --pid
            recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key)
            # Get multiple recipes belonging to project --pid
            for uid in recipes:
                recipe, image = pull_recipe(url=root_url, api_key=api_key, rid=uid)
                write_recipe(recipe=recipe, root_dir=root_dir, image=image)
            return

        project, image = pull_project(pid=pid, url=root_url, api_key=api_key)
        if project:
            write_project(project=project, root_dir=root_dir, image=image)
        return

    def manage_create(self, **options):

        root_dir = options.get("dir")
        json_file = options.get("json")
        job = options.get("create_jobs")
        is_data = options.get("is_data")

        # Get the root from json file.
        if json_file:
            full_path = os.path.abspath(json_file)
            root_dir = root_dir or os.path.dirname(full_path)

        json_files = get_json_files(root_dir=root_dir, json_fname=json_file)

        for fname in json_files:
            json_text = open(os.path.abspath(os.path.join(root_dir, fname)), "r").read()
            recipe_uid = hjson.loads(json_text).get("settings", {}).get("recipe_uid")

            if recipe_uid:
                create_recipe(json_file=fname, root_dir=root_dir, job=job)
            elif is_data:
                data_list = hjson.loads(json_text).get("data", [])
                create_data(root_dir=root_dir, data_list=data_list)
            else:
                create_project(json_file=fname, root_dir=root_dir)

        return

    def manage_list(self, **options):

        url = options.get("url")
        api_key = options.get("key")
        list_ids(url=url, api_key=api_key)
        return

    def add_push_commands(self, parser):

        parser.add_argument('-u', "--url_from_json", action="store_true", help="""Extract url from conf file instead of --url.""")
        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")
        parser.add_argument('--rid', type=str, default="", help="Recipe uid to load.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to load.")
        parser.add_argument('--dir', default='', help="Directory with json files to load in bulk.")
        parser.add_argument('--json', default='', help="""Project or recipe JSON file to load.""")
        return

    def add_pull_commands(self, parser):
        parser.add_argument('-r', "--recipes", action="store_true", help="""Pull recipes of --pid""")
        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")
        parser.add_argument('--rid', type=str, default="", help="Recipe uid to dump.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to dump.")
        parser.add_argument('--dir', default='', help="Directory to store in.")
        return

    def add_list_commands(self, parser):

        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")

        return

    def add_create_commands(self, parser):
        parser.add_argument('--dir', default='', help="Base directory with json files.")
        parser.add_argument('--create_jobs', action="store_true", help="Create job when creating recipe found in --json.")
        parser.add_argument('--is_data', action="store_true", help="Create data object from --json.")
        parser.add_argument('--json', default='', help="""JSON file path relative to create object from.""")
        pass

    def add_arguments(self, parser):

        subparsers = parser.add_subparsers()

        list_parser = subparsers.add_parser("list", help=""" List objects from url or database.""")
        self.add_list_commands(parser=list_parser)

        create_parser = subparsers.add_parser("create", help="""
                                                    Create a project, data or recipe in database from a json file.
                                                    """)
        self.add_create_commands(parser=create_parser)

        push_parser = subparsers.add_parser("push", help="""Update existing project/recipe in remote host or database 
                                                    from a json file.""")
        self.add_push_commands(parser=push_parser)

        pull_parser = subparsers.add_parser("pull", help="""Dump project/recipe into json file from remote 
                                                    host or database.""")
        self.add_pull_commands(parser=pull_parser)

    def handle(self, *args, **options):

        subcommand = sys.argv[2] if len(sys.argv) > 2 else None

        if subcommand == "list":
            self.manage_list(**options)
            return

        if len(sys.argv) <= 3:
            sys.argv.append("--help")
            self.run_from_argv(sys.argv)
            sys.exit()

        if subcommand == "push":
            self.manage_push(**options)
            return
        if subcommand == "create":
            self.manage_create(**options)
            return
        if subcommand == "pull":
            self.manage_pull(**options)
            return
