
from django.core.management.base import BaseCommand


def handle_project():
    return


class Command(BaseCommand):
    help = 'Interact with API end points.'

    def add_arguments(self, parser):
        parser.add_argument('--method', default="GET", help="Request method to use when making request.")
        parser.add_argument('--token', default='', help="API token. Required to access private projects.")
        parser.add_argument("--rid", type=str, default="", help="Recipe uid to dump.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to dump.")
        parser.add_argument("--did", type=str, default="", help="Data uid to dump.")

    def handle(self, *args, **options):
        fname = options['fname']

