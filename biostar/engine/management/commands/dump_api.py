
from django.core.management.base import BaseCommand


def import_from_api(api_url):
    return


class Command(BaseCommand):
    help = 'Import data from api given a base url.'

    def add_arguments(self, parser):

        parser.add_argument('--baseurl', default='', help="Select project by primary id")



        


        pass






