import hjson
import logging
import os

from django.core import management
from django.core.files import File
from django.core.management.base import BaseCommand

from biostar.engine import auth
from biostar.engine.models import Project, User
import sys
from os.path import expanduser

logger = logging.getLogger('engine')


def join(*args):
    return os.path.join(*args)


def parse_json(json, root, privacy=Project.PRIVATE, sticky=False, jobs=False):
    """
    Create a project from a JSON data
    """

    # Load everything relative to the root.
    fpath = os.path.join(root, json)

    # The directory that the project file is located in.
    dirname = os.path.dirname(fpath)

    # This is the
    json_data = hjson.load(open(fpath, 'rb'))

    # The base node of the JSON file.
    base = json_data.get("settings", {})

    if not base:
        logger.error(f"Project {json} must have a 'settings' key")
        sys.exit()

    # The uid is a required element.
    uid = base.get("uid", None)
    if not uid:
        logger.error(f"Project 'settings' dictionary must have a 'uid' field.")
        sys.exit()

    # See if project already exists.
    project = Project.objects.filter(uid=uid).first()

    # Each project uid may only be loaded once.
    if project:
        logger.warning(f"Project uid={project.uid} already exists.")
        return

    # Get more settings into the project.
    name = base.get("name", '')
    text = base.get("text", '')
    summary = base.get("summary", '')
    imgpath = join(dirname, base.get("image", ""))

    # Set the image file stream.
    if os.path.isfile(imgpath):
        stream = File(open(imgpath, 'rb'))
    else:
        stream = None

    # Recipes added at the command line will belong to the superuser.
    user = User.objects.filter(is_superuser=True).first()

    # Setup error. Site has no users.
    if not user:
        logger.error("No valid user was found.")
        return

    # Create the project.
    project = auth.create_project(user=user, uid=uid, summary=summary, name=name, text=text,
                                  stream=stream, privacy=privacy, sticky=sticky)

    # Add extra data specified in the project json file.
    management.call_command("data", json=fpath, id=project.id)

    # Add the analyses specified in the project json.
    analyses = json_data.get("analyses", '')

    # The analyses need are specified relative to the root folder.
    for row in reversed(analyses):
        other_json = os.path.join(root, row['json'])
        template = os.path.join(root, row['template'])
        management.call_command("analysis", id=project.id, add=True, json=other_json, template=template, jobs=jobs)


class Command(BaseCommand):
    help = 'Creates a project.'

    def add_arguments(self, parser):
        parser.add_argument('--root', default='', help="The biostar-recipe project root")
        parser.add_argument('--json', required=True, help="The json file that defines the project relative to the root")
        parser.add_argument('--jobs', action='store_true', default=False, help="Create jobs for the analyses")
        parser.add_argument('--privacy', default="private", help="Privacy of project, defaults to sharable")
        parser.add_argument('--sticky', action='store_true', default=False,
                            help="Make project sticky (high in the order).")

    def handle(self, *args, **options):
        root = options['root']
        json = options['json']
        privacy = options["privacy"].lower()
        sticky = options["sticky"]
        jobs = options["jobs"]

        privacy_map = dict(
            (v.lower(), k) for k, v in Project.PRIVACY_CHOICES
        )

        privacy_value = privacy_map.get(privacy)

        if privacy_value is None:
            logger.error(f"Invalid privacy choice: {privacy}")
            return

        parse_json(json=json, root=root, jobs=jobs, privacy=privacy_value, sticky=sticky)
