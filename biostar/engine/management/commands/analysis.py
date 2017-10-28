import os, sys, logging, hjson, textwrap
from django.core.management.base import BaseCommand
from biostar.engine.models import Job, Project, Analysis, User, Data
from biostar.tools.const import DATA_TYPES
from biostar.engine import auth
from biostar.tools import defaults

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

        # TODO: Impove the help for type
        parser.add_argument('--analysis_type',
                            help=f"Who this job/analysis meant for.",
                            choices=dict(Analysis.TYPE_CHOICES).values())

    def handle(self, *args, **options):

        add = options['add']
        json = options['json']
        pid = options['id']
        template = options['template']
        create_job = options['create_job']

        type_map = lambda dictionary: {y: x for x, y in dictionary.items()}

        analysis_type = type_map(dict(Analysis.TYPE_CHOICES)).get(options.get('analysis_type'), Analysis.USER)

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

            project = Project.objects.get_queryset(user=admin).first()

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
                name = json_data.get("settings", {}).get("name", "No name")
                text = json_data.get("settings", {}).get("help", "No help")
                text = textwrap.dedent(text)
                summary = json_data.get("settings", {}).get("summary", "No summary")
                analysis = auth.create_analysis(project=project,json_text=json_text, summary=summary,
                                                   template=template, name=name, text=text, type=analysis_type)
                logger.info(f"Added analysis '{analysis.name}' to project id={project.id}")

                # Also create a queued job:
                if create_job:
                    # Need to deposit the file as data into the project.
                    # Find all objects that have a path attribute
                    for key, value in json_data.items():
                        path = value.get("path")
                        data_type = value.get("data_type")
                        data_type = DATA_TYPES.get(data_type)

                        if path:
                            data = auth.create_data(project=project,fname=path, data_type=data_type)
                            data.fill_dict(value)
                    auth.create_job(analysis=analysis, json_data=json_data, type=analysis_type)

            except KeyError as exc:
                logger.error(f"processing the analysis: {exc}")
                return
