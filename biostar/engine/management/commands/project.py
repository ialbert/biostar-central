import hjson
import logging
import os

from django.core import management
from django.core.files import File
from django.core.management.base import BaseCommand

from biostar.engine import auth
from biostar.engine.models import Project, User
import sys


logger = logging.getLogger('engine')


class Bunch():
    def __init__(self, **kwargs):
        self.value = ''
        self.name = self.summary = ''
        self.uid = self.text = ''
        self.user = self.stream = None
        self.__dict__.update(kwargs)



def parse_params(fpath, base_dict):
    "Return a Bunch Object of info used to update or create a project."

    # The directory that the project file is located in.
    dirname = os.path.dirname(fpath)

    if not base_dict:
        logger.error(f"Project must have a 'settings' key")
        sys.exit()

    # The uid is a required element.
    uid = base_dict.get("uid", None)
    if not uid:
        logger.error(f"Project 'settings' dictionary must have a 'uid' field.")
        sys.exit()

    imgpath = os.path.join(dirname, base_dict.get("image", ""))

    # Set the image file stream.
    stream = None
    if os.path.isfile(imgpath):
        stream = File(open(imgpath, 'rb'))

    # Recipes added at the command line will belong to the superuser.
    user = User.objects.filter(is_superuser=True).first()

    # Setup error. Site has no users.
    if not user:
        logger.error("No valid user was found.")
        sys.exit()

    bunch = Bunch(uid=uid, name=base_dict.get("name", ''), text= base_dict.get("text", ''),
                  summary= base_dict.get("summary", ''), stream=stream, user=user)
    return bunch


def parse_json(json, root, privacy=Project.PRIVATE, sticky=False, jobs=False, update=False):
    """
    Create a project from a JSON data
    """

    # Load everything relative to the root.
    fpath = os.path.join(root, json)
    json_data = hjson.load(open(fpath, 'rb'))

    # The base node of the JSON file.
    base = json_data.get("settings", {})
    store = parse_params(fpath=fpath, base_dict=base)
    exists = Project.objects.filter(uid=store.uid).exists()
    project = auth.create_project(user=store.user, uid=store.uid, summary=store.summary, name=store.name,
                                text=store.text, stream=store.stream, privacy=privacy, sticky=sticky, update=update)

    # Avoid duplicating data and recipe when updating a project.
    if update or exists:
        action = "updated" if update else "exists"
        logger.info(f"Project {action}. You can update the recipes individually. ")
        return

    # Add extra data specified in the project json file
    management.call_command("data", json=fpath, root=root, id=project.id)

    analyses = json_data.get("analyses", '')
    for row in reversed(analyses):
        json = os.path.join(root, row['json'])
        template = os.path.join(root, row['template'])

        management.call_command("analysis", id=project.id, add=True, json=json, template=template, jobs=jobs)



class Command(BaseCommand):
    help = 'Creates a project.'

    def add_arguments(self, parser):
        parser.add_argument('--root', default='', help="The biostar-recipe project root")
        parser.add_argument('--json', required=True, help="The json file that defines the project relative to the root")
        parser.add_argument('--jobs', action='store_true', default=False, help="Create jobs for the analyses")
        parser.add_argument('--privacy', default="private", help="Privacy of project, defaults to sharable")
        parser.add_argument('--sticky', action='store_true', default=False,
                            help="Make project sticky (high in the order).")
        parser.add_argument('--update', action='store_true', default=False,
                            help="Update an existing project.")

    def handle(self, *args, **options):
        root = options['root']
        json = options['json']
        privacy = options["privacy"].lower()
        sticky = options["sticky"]
        jobs = options["jobs"]
        update = options["update"]

        privacy_map = dict(
            (v.lower(), k) for k, v in Project.PRIVACY_CHOICES
        )

        privacy_value = privacy_map.get(privacy)

        if privacy_value is None:
            logger.error(f"Invalid privacy choice: {privacy}")
            return

        parse_json(json=json, root=root, jobs=jobs, privacy=privacy_value, sticky=sticky, update=update)
