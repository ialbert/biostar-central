import os
from biostar.utils.shortcuts import reverse
from django.conf import settings
from django.core.management.base import BaseCommand


def get_base_url():
    return f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}{settings.HTTP_PORT}"


def export_recipes(base_dir, api_key):

    # Each subdirectory in the base is a project uid

    project_dirs = [os.path.join(base_dir, pid) for pid in os.scandir(base_dir)]

    for fullpath in project_dirs:
        recipes_dir = [os.path.join(fullpath, rid) for rid in os.scandir(fullpath)]

        print(recipes_dir)
    # Each directory in the project uid dir is a recipe



    return


class Command(BaseCommand):
    help = 'Export data from given directory to api using PUT request.'

    def add_arguments(self, parser):

        parser.add_argument('--base', default='', help="Base url to do a reverse look up of api urls.")
        parser.add_argument('--api_key', default='', help="API key used to return all projects.")
        parser.add_argument('--data_dir', default='', help="Directory to export data from.")

        pass

    def handle(self, *args, **options):

        base_url = options["base"] or get_base_url()
        base_dir = options["output"] or settings.API_DUMP
        # Get the project api list
        api_view = reverse("project_api_list")

