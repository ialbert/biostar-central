import logging

import hjson
import os
from django.conf import settings
from django.core.management.base import BaseCommand

from biostar.engine import auth
from biostar.engine.models import Project, Data

logger = logging.getLogger(settings.LOGGER_NAME)


class Bunch():
    def __init__(self, **kwargs):
        self.value = self.uid =  ''
        self.name = self.summary = ''
        self.help = self.type = self.link = ''
        self.__dict__.update(kwargs)


def batch_add_data(data_list, root, update_only, project=None):
    "Iterate over list and add/update data using management commands."

    # Add each collected datatype.
    for bunch in reversed(data_list):
        # This the path to the data.
        path = bunch.value

        # Makes the path relative if necessary.
        path = path if path.startswith("/") else os.path.join(root, path)

        # See if data already exists.
        data = Data.objects.filter(uid=bunch.uid).first()

        if update_only and data:
            # Update existing data
            auth.update_data(data=data, path=path, type=bunch.type,
                             name=bunch.name, summary=bunch.summary, text=bunch.help)
            continue

        if update_only and not data:
            # Skip new data when all we want is to update
            continue

        if (not data) and (not update_only):
            # Create the data if it does not exist.
            auth.create_data(project=project, path=path, type=bunch.type,
                             name=bunch.name, summary=bunch.summary, text=bunch.help)



class Command(BaseCommand):
    help = 'Adds data to a project'

    def add_arguments(self, parser):
        parser.add_argument('--id', default='', help="Select project by primary id")
        parser.add_argument('--uid', default='', help="Select project by unique id")
        parser.add_argument('--json', help="JSON file with the data specification")
        parser.add_argument('--root', default='', help="A root directory for relative paths")
        parser.add_argument('--update_only', default=False, help="Update if uid and update = True in JSON file.")

    def handle(self, *args, **options):

        # Collect the parameters.
        id = options['id']
        uid = options['uid']
        json = options['json']
        root = options['root']
        update_only = options['update_only']

        if not id or uid:
            logger.error(f"Must specify 'id' or 'uid' parameters.")
            return

        # Select project by id or uid.
        if id:
            query = Project.objects.filter(id=id)
        else:
            query = Project.objects.filter(uid=uid)

        # Get the project.
        project = query.first()

        # Project must exist.
        if not project and id:
            logger.error(f"Project with id={id} not found!")
            return

        # Project must exist.
        if not project and id:
            logger.error(f"Project with uid={uid} not found!")
            return

        # There 'data' field of the spec has the files.
        json_data = hjson.load(open(json))
        json_data = json_data.get('data', [])

        # The data field is empty.
        if not json_data:
            logger.error(f"JSON file does not have a valid data field")
            return

        # The datalist is built from the json.
        data_list = [Bunch(**row) for row in json_data]

        # Add each collected datatype.
        batch_add_data(data_list=data_list, root=root, update_only=update_only, project=project)


