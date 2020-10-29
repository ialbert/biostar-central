import hjson
import logging
import sys
import os
from urllib.parse import urljoin
import requests
from django.core.management.base import BaseCommand
from django.shortcuts import reverse

from biostar.recipes.api import tabular_list
from biostar.recipes import auth
from biostar.recipes.models import Project, Analysis

PUSH, PULL = "push", "pull"

logger = logging.getLogger("engine")


def api_response(url, action=PULL, fname=""):
    """
    Send
    """

    # Append token from local settings.TOKEN_FILE.
    params = {"token": auth.get_token()}

    # Push conf as a file to remote host.
    if action == PUSH:
        # Pass conf as a file.
        stream = open(fname, "r")
        files = {"data": stream}
        response = requests.post(url, params=params, files=files)
    # Pull data from remote host.
    else:
        response = requests.get(url, params=params)

    content = response.content.decode()
    return content


def update(obj, data):
    """Resolve the auth function used to update an object."""

    callback = auth.update_project if isinstance(obj, Project) else auth.update_recipe

    return callback(obj=obj, data=data, save=True)


def push_pull(uid, klass, fname="", url=None, action=PULL):
    """
    Push or pull an object from a remote host or database.
    """

    logger.info(f"{action.upper()} action being taken on {url or 'database'}.")

    # Send data to remote host
    if url:
        return api_response(url=url, action=action, fname=fname)

    # Handle object internally.
    obj = klass.objects.filter(uid=uid).first()

    if not obj:
        logger.error("Object does not exist")
        return

    # Push into database
    if action == PUSH:
        stream = open(fname, "r")
        data = hjson.load(stream)
        result = update(obj=obj, data=data)
    # Pull from database
    else:
        result = obj.api_data

    return result


def project_api(uid, fname="", url=None, action=PULL):
    """
    Push or pull a given project from a remote host or database.
    """
    # Construct full API endpoint for project
    endpoint = reverse("project_api", kwargs=dict(uid=uid))

    url = urljoin(url, endpoint) if url else None

    # Push or pull the project
    response = push_pull(klass=Project, uid=uid, fname=fname, url=url, action=action)

    return response


def recipes_api(uid, fname="", url=None, action=PULL):
    """
    Push or pull a given recipe from a remote host or database.
    """

    # Construct full API endpoint for recipe
    endpoint = reverse("recipe_api", kwargs=dict(uid=uid))

    url = urljoin(url, endpoint) if url else None

    # Push or pull the
    response = push_pull(klass=Analysis, uid=uid, fname=fname, url=url, action=action)

    return response


def list_ids(url=None):
    """
    Returns a listing of all project and recipe ids.
    """
    if url:
        endpoint = reverse("api_list")
        url = urljoin(url, endpoint)
        params = {"token": auth.get_token()}
        response = requests.get(url=url, params=params)
        text = response.text
    else:
        text = tabular_list()

    return text


class Command(BaseCommand):
    help = 'Interact with API end points.'

    def add_arguments(self, parser):

        # Give one file or a list of file.
        parser.add_argument('files', metavar='F', default="", nargs='+',
                            help="List of files to push or pull.")

        parser.add_argument('--url', default="",
                            help="URL to get API endpoints from.")

        parser.add_argument("--project", type=str, default="",
                            help="Project uid to push")

        parser.add_argument("--recipe", type=str, default="",
                            help="Recipe uid to push or pull.")

        parser.add_argument("--data", type=str, default="",
                            help="Data uid to push or pull.")

        parser.add_argument("-l", "--list", action="store_true", default=False,
                            help="Show tabular list of projects and recipes..")

        parser.add_argument("-p", "--push", action="store_true", default=False,
                            help="Push data to given project or recipe.")

    def handle(self, *args, **options):
        url = options.get('url')
        pid = options.get('project')
        rid = options.get('recipe')
        fname = options.get('data')
        push = options.get('push')
        files = options.get('push', [])
        listing = options.get('list')

        # List objects
        if listing:
            print(list_ids(url=url))
            return

        if push and not fname:
            logger.error("--data required with --push flag.")
            sys.argv.append("--help")
            self.run_from_argv(sys.argv)
            sys.exit()

        # Resolve what action to is to be taken.
        action = PUSH if push else PULL

        # Call the project API with given pid.
        if pid:
            print(project_api(uid=pid, fname=fname, url=url, action=action))

        # Call the recipe API with given rid.
        if rid:
            print(recipes_api(uid=rid, fname=fname, url=url, action=action))
