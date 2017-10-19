import logging
import os
import sys

import hjson
from django.core.management.base import BaseCommand

from engine.models import Project, User, make_html

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
                json_text = open(spec).read()
                json_data = hjson.loads(json_text)
            except Exception as exc:
                stop(f"error reading out the spec: {exc}")
                json_text = ''
                json_data = {}

            template = open(template).read()

            title = json_data.get("settings", {}).get("title", "No title set")
            text = json_data.get("settings", {}).get("help", "No help set")
            summary = json_data.get("settings", {}).get("summary", "No summary set")
            summary = make_html(summary)
            analysis = project.create_analysis(json_text=json_text, summary=summary,
                                               template=template, title=title, text=text)
            logger.info(f"added analysis {analysis.title} in project {project.title}")
