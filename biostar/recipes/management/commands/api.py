import toml
import sys
import os
from urllib.parse import urljoin
import requests
from django.core.management.base import BaseCommand
from django.shortcuts import reverse

from biostar.recipes.api import tabular_list
from biostar.recipes import auth, models


def api_response(base, endpoint, data={}, method="GET", conf=None):
    """
    Send
    """
    url = urljoin(base, endpoint)

    # Append token from local settings.TOKEN_FILE, if found.
    extra = {"token": auth.get_token()}
    data.update(extra)

    # Push conf as a file to remote host.
    if method == "POST":
        # Pass conf as a file.
        files = {"conf": open(conf, 'r')}
        response = requests.post(url, data=data, files=files)
    # Pull data from remote host.
    else:
        response = requests.get(url, data=data)

    content = response.content.decode()
    return content


def handle_project(uid, conf, url=None, method="GET"):
    """
    Updates an existing project in remote site or local database.
    """

    # Handle the project remotely
    if url:
        endpoint = reverse("project_api", kwargs=dict(uid=uid))
        return api_response(base=url, endpoint=endpoint, method=method, conf=conf)

    # Handle project internally.
    project = models.Project.objects.filter(uid=uid).first()

    # Push conf file into database when passing POST
    if method == "POST":
        data = toml.loads(conf)
        result = auth.update_project(project=project, data=data, save=True)
    # Pull from data base when passing anything else.
    else:
        result = project.api_data

    return result


def handle_recipe(uid, conf, url=None, method="GET"):
    """
    Updates an existing project in remote site or local database.
    """
    # Handle the recipe remotely
    if url:
        endpoint = reverse("recipe_api", kwargs=dict(uid=uid))
        return api_response(base=url, endpoint=endpoint, method=method, conf=conf)

    # Handle recipe internally.
    recipe = models.Analysis.objects.filter(uid=uid).first()

    # Load conf file into database when passing POST
    if method == "POST":
        data = toml.loads(conf)
        result = auth.update_recipe(recipe=recipe, data=data, save=True)
    else:
        result = recipe.api_data

    return result


def list_ids(base=None):
    """
    Returns a listing of all project and recipe ids.
    """
    if base:
        endpoint = reverse("api_list")
        url = urljoin(base, endpoint)
        params = {"token": auth.get_token()}
        response = requests.get(url=url, params=params)
        text = response.text
    else:
        text = tabular_list()

    print(text)


def data_api():
    return


class Command(BaseCommand):
    help = 'Interact with API end points.'

    def list_arguments(self, parser):

        parser.add_argument('--url', default="", help="Site url.")
        return

    def push_arguments(self, parser):

        parser.add_argument('--url', default="", help="URL to get API endpoints from .")
        parser.add_argument("--rid", type=str, default="", help="Recipe uid to push conf file to.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to push conf file to.")
        parser.add_argument("--did", type=str, default="", help="Data uid to push conf file to.")
        parser.add_argument("--conf", type=str, default="",
                            help="Config file with the corresponding project/recipe data to push.")

        return

    def pull_arguments(self, parser):
        parser.add_argument('--url', default="", help="URL to get API endpoints from.")
        parser.add_argument("--rid", type=str, default="", help="Recipe uid to pull.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to pull.")
        parser.add_argument("--did", type=str, default="", help="Data uid to pull.")
        return

    def add_arguments(self, parser):
        subparsers = parser.add_subparsers()

        list_parser = subparsers.add_parser("list", help="""List objects from remote host or database.""")
        push_parser = subparsers.add_parser("push",  help="""Update project/recipe in remote host or database.""")
        pull_parser = subparsers.add_parser("pull", help="""Dump project/recipe from remote host or database.""")
        self.push_arguments(parser=push_parser)
        self.list_arguments(parser=list_parser)
        self.pull_arguments(parser=pull_parser)

    def listing(self, url):
        list_ids(base=url)
        return

    def push_or_pull(self, action, pid, rid, url, conf):

        mapper = {"push": "POST", "pull": "GET"}
        method = mapper.get(action)

        conf = os.path.abspath(conf) if os.path.isfile(conf) else None
        # Handle project API commands
        if pid:
            print(handle_project(url=url, uid=pid, conf=conf, method=method))
            return

        # Handle recipe API commands
        if rid:
            print(handle_recipe(url=url, uid=rid, conf=conf, method=method))
            return

    def handle(self, *args, **options):

        subcommand = sys.argv[2] if len(sys.argv) > 2 else None
        url = options.get('url')
        rid = options.get('rid')
        pid = options.get('pid')
        conf = options.get('conf')

        if len(sys.argv) <= 3:
            sys.argv.append("--help")
            self.run_from_argv(sys.argv)
            sys.exit()

        if subcommand == "list":
            self.listing(url=url)
            return

        self.push_or_pull(action=subcommand, pid=pid, rid=rid, url=url, conf=conf)
