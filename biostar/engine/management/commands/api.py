import logging
import hjson
import os
from urllib.parse import urljoin
from urllib.request import urlopen
import requests
import sys
from functools import partial

from django.conf import settings
from django.core.management.base import BaseCommand
from django.db.models import Q
from django.shortcuts import reverse
from biostar.engine.models import Analysis, Project
from biostar.engine.api import change_image, get_thumbnail
from biostar.engine import auth
from biostar.accounts.models import User

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)


def build_api_url(root_url, uid=None, view="recipe_api_list", api_key=None):

    url = reverse(view, kwargs=dict(uid=uid)) if uid else reverse(view)
    #TODO: use urllib to correctly add params to url
    full_url = urljoin(root_url, url) + f"?k={api_key}"

    return full_url


def remote_upload(stream, root_url, uid, api_key, view):
    """
    Upload data found in stream to root_url.
    Currently uses PUT requests
    """

    payload = dict(k=api_key)

    # Build api url then send PUT request.
    full_url = build_api_url(root_url=root_url, api_key=api_key, view=view, uid=uid)
    response = requests.put(url=full_url, files=dict(file=stream), data=payload)
    if response.status_code == 404:
        print(f"*** Object id : {uid} does not exist on remote host.")
        sys.exit()

    return response


def remote_download(root_url, api_key, view, uid, is_image, outfile, is_json):
    """
    Download data found in root_url using GET request.
    """
    mode = "wb" if is_image else "w"
    # Get data from the api url
    fullurl = build_api_url(root_url=root_url, api_key=api_key, view=view, uid=uid)
    data = urlopen(url=fullurl).read()
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
        name = conf["settings"].get("name", project.name)
        text = conf["settings"].get("text", project.text)
        Project.objects.get_all(uid=uid).update(name=name, text=text)
        return project

    def recipe():
        recipe = Analysis.objects.get_all(uid=uid).first()
        if not recipe:
            # Create empty object if not present then populate.
            project = Project.objects.get_all(uid=pid).first()
            recipe = auth.create_analysis(project=project, json_text="", template="", uid=uid, name="Recipe Name")
        if is_json:
            data = hjson.loads(stream.read())
            name = data["settings"].get("name", recipe.name)
            text = data["settings"].get("help", recipe.text)
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
        recipes = list(filter(filter_func, recipes))
    else:
        query = Q(uid=rid, project__uid=pid) if rid else Q(project__uid=pid)
        recipes = Analysis.objects.get_all().filter(query).values_list("uid", flat=True)

    return recipes


def get_image_name(uid, root_url=None, json="json.hjson", root_dir=None, api_key="", view="recipe_api_json",
                   mtype=Analysis):

    # Get json from url
    if root_url:
        fullurl = build_api_url(root_url=root_url, api_key=api_key, view=view, uid=uid)
        json_text = urlopen(url=fullurl).read().decode()
    # Get json from a file
    elif root_dir:
        path = os.path.join(root_dir, uid, json)
        json_text = open(path).read() if os.path.exists(path) else ""
    # Get json from database
    else:
        json_text = mtype.objects.get_all(uid=uid).first().json_text

    json_settings = hjson.loads(json_text).get("settings", {})

    # Get the local image name from "settings" in json.
    # Defaults to uid.png
    name = json_settings.get("image", f"{uid}.png")

    return name


def recipe_loader(project_dir, pid, api_key="", root_url=None, rid=""):
    """
        Load recipes into api/database from a project found in project_dir.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """
    if not os.path.exists(project_dir):
        print(f"*** Project directory: {project_dir} does not exist.")
        sys.exit()

    # Every subdir in 'project_dir' is a recipe_dir.
    recipe_dirs = [r.name for r in os.scandir(project_dir) if r.is_dir()]
    # Get the specific recipe to load if given.
    recipe_dirs = list(filter(lambda recipe_uid: recipe_uid == rid, recipe_dirs)) if rid else recipe_dirs

    # Prepare the main function used to load.
    load = partial(upload, root_dir=project_dir, root_url=root_url, api_key=api_key, pid=pid, load_recipe=True)

    # Get image name from conf file in directory
    img = lambda uid: get_image_name(uid=uid, root_dir=project_dir)
    for recipe_uid in recipe_dirs:
        load(uid=recipe_uid, fname="json.hjson", view="recipe_api_json", is_json=True)
        load(uid=recipe_uid, fname=img(uid=recipe_uid), view="recipe_api_image", is_image=True)
        load(uid=recipe_uid, fname="template.sh")

        print(f"Loaded recipe id: {recipe_uid}")

    return recipe_dirs


def recipe_dumper(project_dir, pid, root_url=None, api_key="", rid=""):
    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes from API or database.
    recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key, rid=rid)

    dump = partial(download, root_url=root_url, root_dir=project_dir, api_key=api_key)

    # Get image name from json on remote host or local database
    img = lambda uid: get_image_name(uid=uid, root_url=root_url, api_key=api_key)
    for recipe_uid in recipes:
        # Dump json, template, and image for a given recipe
        dump(uid=recipe_uid, fname="json.hjson", is_json=True, view="recipe_api_json")
        dump(uid=recipe_uid, fname=img(uid=recipe_uid), is_image=True, view="recipe_api_image")
        dump(uid=recipe_uid, fname="template.sh")

        print(f"Dumped recipe id: {recipe_uid}")
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

    print(f"Loaded project ({pid}).")
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

    print(f"Dumped project {pid}: {root_dir}.")
    return


class Command(BaseCommand):
    help = 'Dump and load data using api.'

    def add_arguments(self, parser):

        # Load or dump flags
        parser.add_argument('-l', "--load", action="store_true",
                            help="""Load project to url from a directory.
                                    Load to database if --url is not set.""")
        parser.add_argument('-d', "--dump", action="store_true",
                            help="""Dump project from a url to directory. 
                                    Dump from database if --url is not set.""")
        parser.add_argument('--recipes', action="store_true", help="Load/dump recipes from --pid.")
        parser.add_argument('--data', action="store_true", help="""Load/dump data from --pid.
                                                                    Only works in local database.""")

        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")
        parser.add_argument('--dir', default='', help="Directory to store/load data from.")

        parser.add_argument('--pid', type=str, help="Project uid to load or dump.")
        parser.add_argument('--rid', type=str, default="", help="Recipe uid to load or dump.")

    def handle(self, *args, **options):

        load = options.get("load")
        dump = options.get("dump")
        root_url = options["url"]
        api_key = options["key"]
        root_dir = options["dir"] or os.getcwd()
        rid = options["rid"]
        pid = os.path.basename(os.path.dirname(options["pid"]))
        rec = options["recipes"]

        if len(sys.argv) == 2 or not pid:
            self.print_help(sys.argv[0], sys.argv)
            print("\n*** --pid needs to be set.")
            return

        if not (load or dump):
            self.print_help(sys.argv[0], sys.argv)
            print("\n*** Set load (-l) or dump (-d) flag.")
            return

        if load and dump:
            print("\n*** Only one flag can be set.")
            return

        if (root_url and load) and not api_key:
            self.print_help(sys.argv[0], sys.argv)
            print("\n*** --key is required when loading data to remote site.")
            return

        project_dir = os.path.join(root_dir, pid)

        if load:
            # Load the project info and image
            project_loader(pid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key)
            # Load recipes if requested.
            if rec or rid:
                loaded = recipe_loader(project_dir=project_dir, root_url=root_url, api_key=api_key, rid=rid, pid=pid)
                view = reverse("project_view", kwargs=dict(uid=pid))
                view = reverse("recipe_view", kwargs=dict(uid=rid)) if rid else view
                print(f"{len(loaded)} recipes loaded into {urljoin(root_url, view) if root_url else 'database'}")

        elif dump:
            # Dump the project info and image
            project_dumper(pid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key)
            # Dump recipes if requested.
            if rec or rid:
                dumped = recipe_dumper(root_url=root_url, api_key=api_key, project_dir=project_dir, pid=pid, rid=rid)
                print(f"{len(dumped)} recipes dumped into {project_dir}")
