import logging, os, csv
import shutil

from django.conf import settings
from django.core.management.base import BaseCommand
from biostar.engine.models import Data, Job, Analysis, Project

logger = logging.getLogger('engine')

__CURR_DIR = os.path.dirname(os.path.realpath(__file__))


class Command(BaseCommand):
    help = 'Adds access to a project'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):

        data = Data.objects.get_deleted()
        abspath = lambda p: os.path.abspath(p)
        jobs = Job.objects.get_deleted()
        recipes = Analysis.objects.get_deleted()
        projects = Project.objects.get_deleted()

        root_dir = abspath(settings.MEDIA_ROOT)

        def rmdirs(objs):
            for obj in objs:
                obj_path = abspath(obj.get_data_dir())
                # Only delete job files in the media root
                if obj_path.startswith(root_dir):
                    shutil.rmtree(obj_path)
                    print(obj_path, "deleted")
            objs.delete()

        # Delete files associated with objects
        rmdirs(objs=jobs)
        rmdirs(objs=data)

        # Delete files associated with jobs
        recipes.delete()
        projects.delete()

