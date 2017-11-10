import hjson
import logging

from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.engine import auth, const
from biostar.engine.models import Project

logger = logging.getLogger(settings.LOGGER_NAME)


class Bunch():
    def __init__(self, **kwargs):
        self.name = self.path = self.summary = ''
        self.text = self.data_type = self.link = ''
        self.__dict__.update(kwargs)


class Command(BaseCommand):
    help = 'Adds data to a project'

    def add_arguments(self, parser):
        parser.add_argument('--id', default=0, help="Selects project by primary id")
        parser.add_argument('--uid', default="hello", help="Selects project by unique id")
        parser.add_argument('--path', help="The path to the data", default='')
        parser.add_argument('--summary', help="Summary for the data", default='')
        parser.add_argument('--name', help="Name for the data", default='')
        parser.add_argument('--type', help="Data type", default='')
        parser.add_argument('--link', action="store_true", default=False,
                            help="Link the file, not copy")

        parser.add_argument('--json', help="Reads data specification from a json file", default='')

    def handle(self, *args, **options):

        id = options['id']
        uid = options['uid']
        path = options['path'].rstrip("/")
        link = options['link']
        name = options['name']
        summary = options['summary']
        data_type = options['type']

        json = options['json']

        # Select project by id or uid.
        if id:
            query = Project.objects.filter(id=id)
        else:
            query = Project.objects.filter(uid=uid)

        # Get the project.
        project = query.first()

        # Project must exist.
        if not project:
            logger.error(f"Project does not exist: id={id} uid={uid}")
            return

        # Reads a file directly or a spec.
        if not (path or json):
            logger.error(f"Must specify a value for --path or --json")
            return

        if json:
            # There 'data' field of the spec has the files.
            json_data = hjson.load(open(json))
            json_data = json_data.get('data', [])
            data_list = [Bunch(**row) for row in json_data]
        else:
            # There was one data loading request.
            data_list = [
                Bunch(data_type=data_type, path=path, name=name, link=link, summary=summary, text='')
            ]

        for bunch in reversed(data_list):

            # Figure out the bunch datatype
            type_value = const.DATA_TYPE_SYMBOLS.get(bunch.data_type)
            if data_type and not type_value:
                logger.warning(f"Invalid data type: {bunch.data_type}")

            auth.create_data(project=project, path=bunch.path,
                             data_type=type_value,
                             name=bunch.name, link=bunch.link,
                             summary=bunch.summary, text=bunch.text)
