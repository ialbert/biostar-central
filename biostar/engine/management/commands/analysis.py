import os, sys, logging, hjson, textwrap
from django.core.management.base import BaseCommand
from biostar.engine.models import Job, Project, Analysis, User, Data
from biostar.tools import const
from biostar.engine import auth
from biostar.engine import tasks

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

        parser.add_argument('--jobs', action='store_true', default=False,
                            help="Also creates a queued job for the analysis")


    def handle(self, *args, **options):

        add = options['add']
        json = options['json']
        pid = options['id']
        template = options['template']
        jobs = options['jobs']

        admin = User.objects.filter(is_staff=True).first()
        if not admin:
            logger.error("Site has no admin users")
            return

        if not add:
            logger.error("Command requires at least one action: --add --delete")
            return

        if add:

            if not (json and template):
                logger.error("This command requires --json and a --template to be set")
                return

            # Get the target project.
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
                json_path = os.path.dirname(json)
                json_data = hjson.loads(json_text)
            except Exception as exc:
                logger.error(f"Error reading the json: {exc}")
                return

            try:
                # Read the specification
                template = open(template).read()
            except Exception as exc:
                logger.error(f"Error reading template: {exc}")
                return

            try:
                name = json_data.get("settings", {}).get("name", "No name")
                text = json_data.get("settings", {}).get("help", "No help")
                uid = json_data.get("settings", {}).get("uid", "")
                image = json_data.get("settings", {}).get("image", "")
                text = textwrap.dedent(text)
                summary = json_data.get("settings", {}).get("summary", "No summary")

                # Create the analysis
                analysis = auth.create_analysis(project=project, uid=uid, json_text=json_text, summary=summary,
                                                   template=template, name=name, text=text)

                # Load the image if specified.
                if image:
                    image_path = os.path.join(json_path, image)
                    if os.path.isfile(image_path):
                        stream = open(image_path, 'rb')
                        analysis.image.save(image, stream, save=True)
                        logger.info(f"Added image path: {image_path}")
                    else:
                        logger.error(f"Missing image path: {image_path}")

                # Also create a queued job:
                if jobs:
                    # Need to deposit the file as data into the project.
                    # Find all objects that have a path attribute
                    for key, value in json_data.items():
                        path = value.get("path", '')
                        summary = value.get("summary", '')
                        text = value.get("text", '')
                        link = value.get("link", '')
                        data_type = value.get("data_type")
                        name = value.get("name", '') or os.path.basename(path)
                        data_type = const.DATA_TYPE_SYMBOLS.get(data_type)
                        if path:
                            data = auth.create_data(project=project, name=name, path=path, data_type=data_type, link=link,
                                                    summary=summary, text=text)
                            data.fill_dict(value)

                    job_name = f'Results for: {analysis.name}'
                    job = auth.create_job(analysis=analysis, name=job_name, json_data=json_data)

            except KeyError as exc:
                logger.error(f"processing the analysis: {exc}")
                return
