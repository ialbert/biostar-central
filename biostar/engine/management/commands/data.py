import logging, hjson, os
from django.core.management.base import BaseCommand
from biostar.engine.models import Project, User, Data
from biostar.engine import auth
from django.core.files import File
from django.core import management

logger = logging.getLogger('engine')

class Command(BaseCommand):

    help = 'Adds data to a project'

    def add_arguments(self, parser):
        parser.add_argument('--id',  default=0, help="Selects project by primary id")
        parser.add_argument('--uid', default="hello", help="Selects project by unique id")

        parser.add_argument('--path', help="The json file that described the project")

        parser.add_argument('--link', action="store_true", default=False,
                            help="Link the file, not copy")

    def handle(self, *args, **options):
        id = options['id']
        uid = options['uid']

        path = options['path']
        link = options['link']

        if id:
            query = Project.objects.filter(id=id)
        else:
            query = Project.objects.filter(uid=uid)

        project = query.first()

        if not project:
            logger.error(f"Selectd project does not exist: id={id} uid={uid}")
            return
        if not path:
            logger.error(f"Path is required.")
            return
        auth.create_data(project=project, fname=path, link=link)

