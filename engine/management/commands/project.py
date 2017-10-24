from django.core.management.base import BaseCommand
from biostar.tools import defaults
from engine.models import Project

def create():

    pass




class Command(BaseCommand):
    help = 'Creates a project.'

    def add_arguments(self, parser):


        parser.add_argument('--add', action='store_true', default=False,
                            help="Adds an analysis to the project")

        parser.add_argument('--json',
                            help="The json specification file ( if ")
        parser.add_argument('--template',
                            help="The template for the analysis")

        parser.add_argument('--create_job', action='store_true', default=False,
                            help="Also creates a queued job for the analysis")
        # TODO: Impove the help for usage

        parser.add_argument('--project_usage',
                            help=f"Who this job/analysis meant for.",
                            default=defaults.USAGE, choices=dict(Project.USAGE_CHOICES).values())

    def handle(self, *args, **options):
        add = options['add']
        json = options['json']
        pid = options['id']
        template = options['template']
        create_job = options['create_job']

        usage_map = lambda dictionary: {y: x for x, y in dictionary.items()}
        project_usage = usage_map(dict(Project.USAGE_CHOICES)).get(options['project_usage'], Project.USER)
        print(project_usage, "SSSSS")

        Project.objects.filter(id=pid).update(usage=project_usage)
