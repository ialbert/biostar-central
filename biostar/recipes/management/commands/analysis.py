import logging
import os
import textwrap

import toml as hjson
from django.core.management.base import BaseCommand

from biostar.recipes import auth
from biostar.recipes.models import Project, Analysis, Data

logger = logging.getLogger('engine')

__CURR_DIR = os.path.dirname(os.path.realpath(__file__))


class Command(BaseCommand):
    help = 'Manages analyses.'

    def add_arguments(self, parser):

        parser.add_argument('--add', action='store_true', default=False,
                            help="Adds an analysis to a project")

        parser.add_argument('--id',
                            help="Specifies the project id")

        parser.add_argument('--uid',
                            help="Specifies the project uid")

        parser.add_argument('--json',
                            help="The json specification file")

        parser.add_argument('--template',
                            help="The template for the analysis")

        parser.add_argument('--jobs', action='store_true', default=False,
                            help="Also creates a queued job for the analysis")

        parser.add_argument('--update', action='store_true', default=False,
                            help="Update an existing recipe")

    def handle(self, *args, **options):

        json = options['json']
        pid = options['id']
        uid = options["uid"]
        template_fname = options['template']
        jobs = options['jobs']
        update = options["update"]

        verbosity = int(options['verbosity'])

        if verbosity > 1:
            logger.setLevel(logging.DEBUG)
            logger.info(f"level={verbosity}")

        # Require JSON and templates to exist.
        if not (json and template_fname):
            logger.error("This command requires --json and a --template to be set")
            return

        # Get the target project.
        if pid:
            project = Project.objects.filter(id=pid).first()
        else:
            project = Project.objects.filter(uid=uid).first()

        # Invalid project specified.
        if not project:
            logger.error(f'No project with id={pid} , uid={uid}')
            return

        # JSON file does not exist.
        if not os.path.isfile(json):
            logger.error(f'No file found for --json={json}')
            return

        # Template file does not exist.
        if not os.path.isfile(template_fname):
            logger.error(f'No file found for --template={template_fname}')
            return

        try:
            # Parse the json_text into json_data
            json_text = open(json).read()
            json_path = os.path.dirname(json)
            json_data = hjson.loads(json_text)
        except Exception as exc:
            logger.exception(f"JSON exception in file: {json}\n{exc}")
            return
        try:
            # Read the specification
            template = open(template_fname).read()
        except Exception as exc:
            logger.exception(f"Template exception: {exc}")
            return

        try:
            name = json_data.get("settings", {}).get("name", "")
            text = json_data.get("settings", {}).get("help", "")
            uid = json_data.get("settings", {}).get("uid", "")
            image = json_data.get("settings", {}).get("image", "")
            text = textwrap.dedent(text)
            summary = json_data.get("settings", {}).get("summary", "")

            analysis = auth.create_analysis(
                project=project, uid=uid, json_text=json_text,
                summary=summary, template=template, name=name, text=text, security=Analysis.AUTHORIZED, update=update)

            # Load the image if specified.
            if image:
                image_path = os.path.join(json_path, image)
                if os.path.isfile(image_path):
                    stream = open(image_path, 'rb')
                    analysis.image.save(image, stream, save=True)
                    logger.info(f"Image path: {image_path}")
                else:
                    logger.error(f"Skipping invalid image path: {image_path}")

            if jobs:
                # When creating a job automatically for data in projects
                # it will try to match the value of the parameter to the data name.
                missing_name = ''
                for key, obj in json_data.items():
                    if obj.get("source") != "PROJECT":
                        continue

                    name = obj.get('value', '')
                    data = Data.objects.filter(project=project, name=name).first()

                    if not data:
                        missing_name = name
                        break

                    data.fill_dict(obj)

                if missing_name:
                    logger.error(f"Job not created! Missing data:{missing_name} in analysis:{analysis.name}")
                else:
                    auth.create_job(analysis=analysis, json_data=json_data)

        except Exception as exc:
            logger.exception(f"Error: {exc}")
            return
