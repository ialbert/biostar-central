import hjson
import logging
import os

from django.core import management
from django.core.files import File
from django.core.management.base import BaseCommand

from biostar.engine import auth
from biostar.engine.models import Project, User

logger = logging.getLogger('engine')


def join(*args):
    return os.path.join(*args)


def parse_json(json, privacy=Project.PRIVATE, sticky=False, jobs=False):
    """
    Create a project from a JSON data
    """
    data = hjson.load(open(json, 'rb'))
    dirname = os.path.dirname(json)

    email = data.get("email", 'email')
    uid = data.get("uid", None)
    name = data.get("name", '')
    text = data.get("text", '')
    summary = data.get("summary", '')
    imgpath = join(dirname, data.get("image", ""))

    project = Project.objects.filter(uid=uid).first()

    if project:
        logger.info(f"Project uid={project.uid} already exists")
        return

    if os.path.isfile(imgpath):
        stream = File(open(imgpath, 'rb'))
    else:
        stream = None

    user = User.objects.filter(email=email).first()
    user = user or User.objects.filter(is_superuser=True).first()

    project = auth.create_project(user=user, uid=uid, summary=summary, name=name, text=text,
                                  stream=stream, privacy=privacy, sticky=sticky)

    analyses = data.get("analyses", '')

    for row in analyses:
        ajson = row['json']
        template = row['template']
        management.call_command("analysis", id=project.id, add=True, json=ajson, template=template, create_job=jobs)

    print (json)
    management.call_command("data", json=json, id=project.id)


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
