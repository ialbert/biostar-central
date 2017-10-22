from django.core.management.base import BaseCommand
from engine.models import Data, Project, get_datatype


def copy_data(data_id, project_id=1, fname=None):


    data = Data.objects.filter(id=data_id).first()


    project = Project.object.filter(id=project_id).first()


    project.create_data(stream=stream, name=name, data_type=data_type, text=text,
                        owner=owner)

    # get the name from the file name
    return




class Command(BaseCommand):
    help = 'Upload Data'

    def add_arguments(self, parser):

        parser.add_argument('--data_id', required =True,
                            help="Upload file to data list of project(s).")
        parser.add_argument('--ids', default=1,
                            help="One or more project ids (comma separated)")
        parser.add_argument('--path', default=False, action="store_true",
                            help="New filepath of uploaded data.")
        parser.add_argument('--path', default=False, action="store_true",
                            help="New filepath of uploaded data.")


    def handle(self, *args, **options):

            fname = options["fname"]
            ids = options['ids']
            data_path = options['path']
