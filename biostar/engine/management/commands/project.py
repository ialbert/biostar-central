import logging
from django.core.management.base import BaseCommand
from biostar.tools import defaults
from biostar.engine.models import Project, Analysis, Data, User
from django.core import management

logger = logging.getLogger('engine')

def create(owner, name=defaults.DEMO_PROJECT_NAME, summary=defaults.DEMO_PROJECT_SUMMARY,
           add=False, json=None, template=None, create_job=False, project_type=None, analysis_type=None):

    if not owner.is_superuser:
        logger.error(f'User is not admin.')
        return

    project = Project.objects.create(owner=owner, name=name, summary=summary, type=project_type)
    logger.info(f'Created project name={project.name}, id={project.id} with type:{project.get_type_display()}')
    project.save()

    if add:
        assert json and template
        assert isinstance(json, str) and isinstance(template, str)
        management.call_command('analysis', template=template, id=project.id, create_job=create_job, json=json,
                                type=analysis_type, add=True)


class Command(BaseCommand):

    help = 'Creates a project.'

    def add_arguments(self, parser):

        parser.add_argument('--name',
                            help=f"Name of created project ( in quotes). default = {defaults.PROJECT_NAME}",
                            default=defaults.PROJECT_NAME)

        parser.add_argument('--summary',
                            help=f"Name of created project ( in quotes). default = {defaults.PROJECT_SUMMARY}",
                            default=defaults.PROJECT_SUMMARY)

        parser.add_argument('--creator_email',
                            help=f"Name of created project. default = first admin user")

        parser.add_argument('--add', action='store_true', default=False,
                            help="Adds an analysis to the project")

        parser.add_argument('--json',
                            help="The json specification file ( if ")

        parser.add_argument('--template',
                            help="The template for the analysis")

        parser.add_argument('--create_job', action='store_true', default=False,
                            help="Also creates a queued job for the analysis")
        # TODO: Impove the help for type

        parser.add_argument('--project_type',
                            help=f"Who this job/analysis meant for.",
                            default=defaults.USAGE, choices=dict(Project.TYPE_CHOICES).values())

        parser.add_argument('--analysis_type',
                            help=f"Who this job/analysis meant for.",
                            default=defaults.USAGE, choices=dict(Analysis.TYPE_CHOICES).values())


    def handle(self, *args, **options):

        name = options['name']
        summary = options['summary']

        add = options['add']
        json = options.get('json')
        template = options.get('template')
        create_job = options['create_job']
        type = options['project_type']
        creator_email = options.get('creator_email', '')

        owner = User.objects.filter(is_superuser=True, email=creator_email).first()
        if not owner:
            owner = User.objects.filter(is_superuser=True).first()
        type_map = lambda dictionary: {y: x for x, y in dictionary.items()}
        project_type = type_map(dict(Project.TYPE_CHOICES)).get(type, Project.USER)

        if add:
            assert json and template
        create(owner, name, summary, add, json, template, create_job, project_type)
