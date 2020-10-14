from urllib.parse import urljoin
import requests
from django.core.management.base import BaseCommand
from django.shortcuts import reverse

from biostar.recipes.api import tabular_list


def api_request(url, method="GET", params={}, data={}):

    if method == "POST":
        response = requests.post(url, params=params, data=data)
    else:
        response = requests.get(url, params=params, data=data)

    return response


def project_api(url, uid, token=None, method="GET", data=dict()):
    """
    Update a project at a given url
    """

    endpoint = reverse("project_api", kwargs=dict(uid=uid))
    full_url = urljoin(url, endpoint)
    params = {"token": token}

    response = api_request(url=full_url, method=method, params=params)

    print(url, response)

    return


def recipe_api(url, uid, token=None, method="GET", data=dict()):


    endpoint = reverse("recipe_api", kwargs=dict(uid=uid))
    params = {"token": token}
    full_url = urljoin(url, endpoint)

    response = api_request(url=full_url, method=method, params=params)

    return


def list_ids(url=None, token=""):
    """
    Returns a listing of all project and recipe ids.
    """
    if url:
        endpoint = reverse("api_list")
        response = requests.get(url=endpoint)
        text = response.text
    else:
        text = tabular_list()

    print(text)


def data_api():
    return


class Command(BaseCommand):
    help = 'Interact with API end points.'

    def add_arguments(self, parser):
        parser.add_argument('--method', default="get", help="Request method to use when making request.")
        parser.add_argument('--url', default="GET", help="URL to get API endpoints from.")
        parser.add_argument('--token', default='', help="API token. Required to access private projects.")
        parser.add_argument("--list", type=str, default="", help="List projects and recipes.")
        parser.add_argument("--rid", type=str, default="", help="Recipe uid to dump.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to dump.")
        parser.add_argument("--did", type=str, default="", help="Data uid to dump.")

    def handle(self, *args, **options):
        fname = options['fname']

        url = options['url']
        token = options['token']
        list = options['list']
        method = options['method']
        rid = options['rid']
        pid = options['pid']
        did = options['did']

        if list:
            list_ids(url=url)

        if rid:
            recipe_api(url=url, uid=rid, token=token, method=method)

        if pid:
            project_api(url=url, uid=rid, token=token, method=method)

        #if did:
