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
        pid = options['pid']
        template = options['template']

        admin = User.objects.filter(is_staff=True).first()
        if not admin:
            stop("site has no admin users")

        if not add or 0:
            stop("command requires at least one action: --add --delete")
        if add:

            if not (pid and spec and template):
                stop("the command requires --pid --spec --template to be set")

            project = Project.objects.filter(id=pid).first()
            if not project:
                stop(f'No project with id={pid}')

            if not os.path.isfile(spec):
                stop(f'No file found for --spec={spec}')

            if not os.path.isfile(template):
                stop(f'No file found for --template={template}')

            try:
                json_data = open(spec).read()
                json_data = hjson.loads(json_data)
            except Exception as exc:
                stop(f"error reading out the spec: {exc}")

            template = open(template).read()

            title = json_data.get("settings", {}).get("title", "No title set")
            text = json_data.get("settings", {}).get("help", "No help set")
            analysis = Analysis(owner=admin, project=project, title=title, text=text,
                        json_data=json_data, makefile_template=template)
            analysis.save()
            logger.info(f"added analysis {analysis.title} in project {project.title}")


