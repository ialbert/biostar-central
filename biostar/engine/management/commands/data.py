import logging, os
from django.core.management.base import BaseCommand
from biostar.engine.models import Project
from biostar.engine import auth
from django.conf import settings

logger = logging.getLogger(settings.LOGGER_NAME)

class Command(BaseCommand):

    help = 'Adds data to a project'

    def add_arguments(self, parser):
        parser.add_argument('--id',  default=0, help="Selects project by primary id")
        parser.add_argument('--uid', default="hello", help="Selects project by unique id")
        parser.add_argument('--path', help="The path to the data")
        parser.add_argument('--summary', help="Summary for the data", default='')
        parser.add_argument('--link', action="store_true", default=False,
                            help="Link the file, not copy")

    def handle(self, *args, **options):

        id = options['id']
        uid = options['uid']
        path = options['path'].rstrip("/")
        link = options['link']
        summary = options['summary']

        if id:
            query = Project.objects.filter(id=id)
        else:
            query = Project.objects.filter(uid=uid)

        # Get the project.
        project = query.first()

        if not project:
            logger.error(f"Project does not exist: id={id} uid={uid}")
            return

        if not path:
            logger.error(f"Must specify a value for --path")
            return

        auth.create_data(project=project,  path=path, link=link, summary=summary)
