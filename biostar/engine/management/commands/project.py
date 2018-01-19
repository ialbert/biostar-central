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


def join(*args):
    return os.path.join(*args)


def parse_json(json, privacy=Project.PRIVATE, sticky=False, jobs=False):
    """
    Create a project from a JSON data
    """
    data = hjson.load(open(json, 'rb'))
    dirname = os.path.dirname(json)

    root = data.get("settings", {})

    if not root:
        logger.error(f"Project {json} must have a 'settings' key")
        sys.exit()

    uid = root.get("uid", None)
    name = root.get("name", '')
    text = root.get("text", '')
    summary = root.get("summary", '')
    imgpath = join(dirname, root.get("image", ""))

    project = Project.objects.filter(uid=uid).first()

    # Stop if the project already exists.
    if project:
        logger.info(f"Project uid={project.uid} already exists")
        return

    # Set the image file stream.
    if os.path.isfile(imgpath):
        stream = File(open(imgpath, 'rb'))
    else:
        stream = None

    # Recipes added at the command line belong to the superuser.
    user = User.objects.filter(is_superuser=True).first()

    # Setup error. Site has no users.
    if not user:
        logger.error("No valid users found")
        return

    # Create the project.
    project = auth.create_project(user=user, uid=uid, summary=summary, name=name, text=text,
                                  stream=stream, privacy=privacy, sticky=sticky)

    # Add datatypes specific to a project
    datatypes = data.get("datatypes", '')

    for name in datatypes:
        symbol = datatypes[name].get("symbol", '')
        help = datatypes[name].get("help", '')
        datatype = auth.create_datatype(name=name, symbol=symbol, help=help, project=project)
        datatype.save()

    # Add extra data specified in the project json file.
    management.call_command("data", json=json, id=project.id)

    # Add the analyses specified in the project json.
    analyses = data.get("analyses", '')

    for row in reversed(analyses):
        other_json = row['json']
        template = row['template']
        management.call_command("analysis", id=project.id, add=True, json=other_json, template=template, jobs=jobs)




class Command(BaseCommand):
    help = 'Creates a project.'

    def add_arguments(self, parser):
        parser.add_argument('--json', required=True, help="The json file that described the project")
        parser.add_argument('--jobs', action='store_true', default=False, help="Create jobs for the analyses")
        parser.add_argument('--privacy', default="private", help="Privacy of project, defaults to sharable")
        parser.add_argument('--sticky', action='store_true', default=False,
                            help="Make project sticky (high in the order).")

    def handle(self, *args, **options):
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

        parse_json(json, jobs=jobs, privacy=privacy_value, sticky=sticky)
