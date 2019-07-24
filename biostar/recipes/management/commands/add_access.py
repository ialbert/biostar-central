import logging, os, csv

from django.core.management.base import BaseCommand

from biostar.recipes import models
from biostar.recipes.models import Access

logger = logging.getLogger('engine')

__CURR_DIR = os.path.dirname(os.path.realpath(__file__))

# The valid access choices.
CHOICE_MAP = dict(
    read=Access.READ_ACCESS,
    write=Access.WRITE_ACCESS,
    owner=Access.WRITE_ACCESS
)

CHOICES = list(CHOICE_MAP.keys())


class Command(BaseCommand):
    help = 'Adds access to a project'

    def add_arguments(self, parser):
        parser.add_argument('fname', default=0,
                            help="Specifies project by id")

    def handle(self, *args, **options):
        fname = options['fname']

        if not os.path.isfile(fname):
            logger.error(f'Not a valid filename: {fname}')
            return

        stream = open(fname, 'rU')

        for row in csv.DictReader(stream):
            email = row.get("Email", '')
            uid = row.get("Uid", '')
            access = row.get("Access", '')

            # All fields must be filled in.
            if not (uid and email and access):
                logger.error(f"Invalid row: {row}")
                continue

            # Fix capitalization and white spaces.
            email = email.strip().lower()
            uid = uid.strip()
            access = access.strip().lower()

            # Get the project by UID.
            project = models.Project.objects.filter(uid=uid).first()

            if not project:
                logger.error(f"Project not found: uid={uid}")
                continue

            # Get the user by email.
            user = models.User.objects.filter(email=email).first()

            if not user:
                logger.error(f"User not found: email={email}")
                return

            # Drop all other permissions for the user.
            models.Access.objects.filter(user=user, project=project).delete()

            # Get the access value and create the access entry.
            access_value = CHOICE_MAP.get(access)
            if not access_value:
                logger.error(f"Invalid access type: uid={project.uid} email={email}  access={access}")
                continue

            if access == "owner":
                project.owner = user
                project.save()
            else:
                # Create the access types
                models.Access.objects.create(user=user, project=project, access=access_value)

            logger.info(f'Access created uid={project.uid} email={email} access={access}')
