from django.core.management.base import BaseCommand
from biostar.engine.search import write_to_index
from biostar.engine.models import Project, Analysis, Data, Job


class Command(BaseCommand):
    help = 'Add a batch of objects to the search index'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):

        for modelname in [Project, Analysis, Data, Job]:

            queryset = modelname.objects.get_all()

            # Write objects from given queryset into index
            print(write_to_index(queryset=queryset))




