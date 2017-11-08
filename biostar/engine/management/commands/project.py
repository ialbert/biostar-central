import logging, hjson, os
from django.core.management.base import BaseCommand
from biostar.engine.models import Project, User
from biostar.engine import auth
from django.core.files import File
from django.core import management

logger = logging.getLogger('engine')

def join(*args):
    return os.path.join(*args)

def parse_json(path, privacy=Project.SHAREABLE):
    """
    Create a project from a JSON data
    """
    data = hjson.load(open(path, 'rb'))
    dirname = os.path.dirname(path)

    email = data.get("email", 'email')
    uid = data.get("uid",None)
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
                                  stream=stream, privacy=privacy)

    analyses = data.get("analyses", '')

    for row in analyses:
        json = row['json']
        template = row['template']
        management.call_command("analysis", id=project.id, add=True, json=json, template=template, create_job=True)


class Command(BaseCommand):

    help = 'Creates a project.'

    def add_arguments(self, parser):
        parser.add_argument('--json', required=True, help="The json file that described the project")
        parser.add_argument('--privacy', default="share", help="The json file that described the project")


    def handle(self, *args, **options):
        path = options['json']
        privacy = options["privacy"]


        parse_json(path)






