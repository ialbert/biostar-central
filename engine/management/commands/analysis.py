import os, sys, logging, hjson
from engine.models import Project, User, make_html
from django.core.management.base import BaseCommand

logger = logging.getLogger('engine')

__CURR_DIR = os.path.dirname(os.path.realpath(__file__))

class Command(BaseCommand):
    help = 'Manages analyses.'

    def add_arguments(self, parser):

        parser.add_argument('--add', action='store_true', default=False,
                            help="Adds an analysis to a project")
        parser.add_argument('--id', default=1,
                            help="Specifies the project id")
        parser.add_argument('--json',
                            help="The json specification file")
        parser.add_argument('--template',
                            help="The template for the analysis")
        parser.add_argument('--create_job', action='store_true', default=False,
                            help="Also creates a queued job for the analysis")

    def handle(self, *args, **options):

        add = options['add']
        json = options['json']
        pid = options['id']
        template = options['template']
        create_job = options['create_job']

        admin = User.objects.filter(is_staff=True).first()
        if not admin:
            logger.error("site has no admin users")
            return

        if not add:
            logger.error("command requires at least one action: --add --delete")
            return

        if add:

            if not (json and template):
                logger.error("this command requires --json --template to be set")
                return

            project = Project.objects.filter(id=pid).first()
            if not project:
                logger.error(f'No project with id={pid}')
                return

            if not os.path.isfile(json):
                logger.error(f'No file found for --json={json}')
                return

            if not os.path.isfile(template):
                logger.error(f'No file found for --template={template}')
                return

            try:
                # Parse the json_text into json_data
                json_text = open(json).read()
                json_data = hjson.loads(json_text)
            except Exception as exc:
                logger.error(f"error leading the template: {exc}")
                return

            try:
                # Read the specification
                template = open(template).read()
            except Exception as exc:
                logger.error(f"error reading out the spec: {exc}")
                return

            try:
                title = json_data.get("settings", {}).get("title", "No title set")
                text = json_data.get("settings", {}).get("help", "No help set")
                summary = json_data.get("settings", {}).get("summary", "No summary set")
                summary = make_html(summary)
                analysis = project.create_analysis(json_text=json_text, summary=summary,
                                                   template=template, title=title, text=text)
                logger.info(f"Added analysis '{analysis.title}' to project id={project.id}")

                # Also create a queued job:
                if create_job:
                    # Need to deposit the file as data into the project.
                    # Find all objects that have a path attribute
                    for key, value in json_data.items():
                        path = value.get("path")
                        if path:
                            data = project.create_data(fname=path)
                            json_data[key]['path'] = data.get_path()
                    analysis.create_job(json_data=json_data)

            except KeyError as exc:
                logger.error(f"processing the analysis: {exc}")
                return
