
import logging, uuid
from engine.const import *
import os
from django.conf import settings
from django.core import management
import hjson as json
from django.core.files import File
from . import util


logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


def get_uuid(limit=None):
    return str(uuid.uuid4())[:limit]


def init_proj(sender, **kwargs):
    """
    Populate initial projects with N number data
    Creates one analysis model to allow for jobs to be run
    """
    from engine.models import User, Project, Data, Analysis, Job, Group

    owner = User.objects.all().first()

    # Give Lamar group but project inherits the group from owner so
    # this is overridden
    group = owner.groups.first()

    # Needs to run only if there are no projects.
    if Project.objects.filter().all():
        return

    # Make the test projects.
    from biostar.tools.testdata import TEST_PROJECTS, TEST_DATA, TEST_ANALYSES

    for title, description in TEST_PROJECTS:

        project = Project(title=title, owner=owner, text=description,
                          group=group)
        project.save()

        logger.info(f'Created project: {project.title} belonging to {group} group.')

        # Add data to each project.
        for data_title, data_desc, data_file, data_type in TEST_DATA:
            # data naming is a bit weired.
            stream = File(open(data_file, 'rb'))

            data_file = os.path.split(data_file)[-1]
            size = f"{stream.size}"
            data = Data(title=data_title, owner=owner, text=data_desc, project=project, data_type=data_type, size=size)
            data.file.save(data_file, stream, save=True)
            data.save()

        # Initialize the same analyses for each project.
        for spec_path, tmpl_path in TEST_ANALYSES:
            management.call_command("analysis", add=True, pid=project.id, spec=spec_path, template=tmpl_path)

        # Get the fastqc analysis.
        analysis = Analysis.objects.filter(title__startswith="Generate FastQC").first()

        # Get a FASTQ data
        fastq_data = Data.objects.filter(data_type=FASTQ_TYPE).first()

        # This work for this analysis only.
        json_data = json.loads(analysis.json_text)
        json_data['data']['path'] = fastq_data.get_path()
        json_text = json.dumps(json_data)

        # Create four jobs for each project.
        for state in [Job.ERROR, Job.QUEUED]:
            analysis.create_job(state=state, json_text=json_text)


def init_users(sender, **kwargs):
    """
    Creates admin users if these are not present.
    """
    from engine.models import User, Group

    groups = Group.objects.all()
    logger.info("Setting up users")

    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            user = User(first_name=name, email=email,
                        is_superuser=True, is_staff=True)
            user.set_password(settings.SECRET_KEY)
            user.save()
            logger.info(f"created admin user: {user.email}")

            # add admin to all groups ( for now atleast )
            for group in groups:

                group.user_set.add(user)
                group.save()
                logger.info(f"adding {user.email} to {group} group.")

    for email, user_groups in settings.REGULAR_TEST_USERS.items():

        if not User.objects.filter(email=email):

            test_user = User(email=email, username=util.get_uuid())
            test_user.set_password(email)
            logger.info(f"creating user: {test_user.email}")
            test_user.save()

            for gr in user_groups:
                user_group = Group.objects.filter(name=gr).first()
                test_user.groups.add(user_group)
                user_group.save()
                logger.info(f"adding {test_user.email} to {user_group} group.")



def init_groups(sender, **kwargs):

    from engine.models import Group

    logger.info("Setting up Groups")

    if Group.objects.filter().all():
        return

    for name in INITIAL_GROUPS:

        group = Group.objects.create(name=name)
        group.save()

    return


def init_site(sender, **kwargs):
    """
    Updates site domain and name.
    """
    from django.contrib.sites.models import Site

    # Print information on the database.
    db = settings.DATABASES['default']
    logger.info("db.engine={}, db.name={}".format(db['ENGINE'], db['NAME']))
    logger.info("email.backend={}".format(settings.EMAIL_BACKEND))
    logger.info("default.email={}".format(settings.DEFAULT_FROM_EMAIL))

    # Create the default site if necessary.
    Site.objects.get_or_create(id=settings.SITE_ID)

    # Update the default site domain and name.
    Site.objects.filter(id=settings.SITE_ID).update(domain=settings.SITE_DOMAIN, name=settings.SITE_NAME)

    # Get the current site
    site = Site.objects.get(id=settings.SITE_ID)
    logger.info("site.name={}, site.domain={}".format(site.name, site.domain))


