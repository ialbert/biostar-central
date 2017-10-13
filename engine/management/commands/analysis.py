from django.core.management.base import BaseCommand
from engine.models import Analysis, Project, User
import os, sys,hjson, logging

logger = logging.getLogger('engine')

__CURR_DIR = os.path.dirname(os.path.realpath(__file__))

def stop(msg):
    logger.error(msg)
    sys.exit()

class Command(BaseCommand):
    help = 'Manages analyses.'

    def add_arguments(self, parser):

        # Named (optional) arguments
        parser.add_argument('--add',
                            action='store_true',
                            default=False,
                            required=False, help="Adds an analysis to a project")
        parser.add_argument('--pid',
                            required=False, help="Specifies the project id")
        parser.add_argument('--spec',
                            required=False, help="The hsjon specification file")
        parser.add_argument('--template',
                            required=False, help="The template for the analysis")

    def handle(self, *args, **options):

        add = options['add']
        spec = options['spec']
        project_id = options['pid']
        template = options['template']

        admin = User.objects.filter(is_staff=True).first()
        if not admin:
            stop("site has no admin users")

        if add:

            if not (project_id and spec and template):
                stop("the command requires --project_id --spec --template to be set")

            project = Project.objects.filter(id=project_id).first()
            if not project:
                stop(f'No project with id={project_id}')

            if not os.path.isfile(spec):
                stop(f'No file found for --spec={spec}')

            if not os.path.isfile(template):
                stop(f'No file found for --template={template}')

            json_data = open(spec).read()
            makefile_template = open(template).read()

            analysis = Analysis(owner=admin, project=project,
                        json_data=json_data, makefile_template=makefile_template)
            analysis.save()
            logger.info(f"added analysis {analysis.title} to project {project.title}")


