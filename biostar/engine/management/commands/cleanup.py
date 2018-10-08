import logging, os, csv
import shutil

from django.core.management.base import BaseCommand
from biostar.engine.models import Data, Job

logger = logging.getLogger('engine')

__CURR_DIR = os.path.dirname(os.path.realpath(__file__))


def delete(queryset):

    for obj in queryset:
        basedir = obj.get_data_dir()
        files = os.listdir(basedir)

        files = list(map(os.path.join, [basedir] * len(files), files))
        linked_files = list(filter(lambda p: os.path.islink(p), files))

        # Delete the object
        obj.delete()

        # Leave linked files alone
        if linked_files:
            continue

        # Delete the root dir of the object
        shutil.rmtree(basedir)

        print(basedir, "deleted")


class Command(BaseCommand):
    help = 'Adds access to a project'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):

        data = Data.objects.get_deleted()
        jobs = Job.objects.get_deleted()

        job_paths = list(jobs.values_list('path', flat=True))

        # Delete jobs
        jobs.delete()

        # Delete files associated with jobs
        list(map(shutil.rmtree, job_paths))

        # Delete the data next
        delete(data)
