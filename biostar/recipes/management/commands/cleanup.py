import logging, os, csv
import shutil

from django.conf import settings
from django.core.management.base import BaseCommand
from biostar.recipes.models import Data, Job, Analysis, Project

logger = logging.getLogger('engine')

__CURR_DIR = os.path.dirname(os.path.realpath(__file__))


class Command(BaseCommand):
    help = 'Adds access to a project'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):

        data = Data.objects.filter(deleted=True)
        jobs = Job.objects.filter(deleted=True)
        recipes = Analysis.objects.filter(deleted=True)
        projects = Project.objects.filter(deleted=True)

        root_dir = os.path.abspath(settings.MEDIA_ROOT)

        def rmdirs(objs):
            for obj in objs:
                obj_path = os.path.abspath(obj.get_data_dir())
                # Only delete job files in the media root
                if obj_path.startswith(root_dir):
                    shutil.rmtree(obj_path)
                    logger.info(f"{obj_path} deleted")
            objs.delete()

        # Delete files associated with objects
        rmdirs(objs=jobs)
        rmdirs(objs=data)

        recipes.delete()
        projects.delete()

