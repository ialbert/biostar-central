import logging, uuid
import os
from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate
from engine.const import *
from django.core import management
from copy import copy
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
    from engine.models import Project, Data, Analysis, make_job, Job
    from engine.models import User

    owner = User.objects.all().first()

    # Needs to run only if there are no projects.
    if Project.objects.filter().all():
        return

    # Make the test projects.
    for title, description in TEST_PROJECTS:

        project = Project(title=title, owner=owner, text=description)
        project.save()

        logger.info(f'Created project: {project.title}')

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
            management.call_command("analysis", add=True, pid=project.id,spec=spec_path, template=tmpl_path)

        # Get the fastqc analysis.
        fastq_analysis = Analysis.objects.filter(title__startswith="Fastqc report").first()

        # Get a FASTQ data
        fastq_data = Data.objects.filter(data_type=FASTQ_TYPE).first()

        filled_json = json.loads(fastq_analysis.json_text)
        filled_json['data']['path'] = fastq_data.file.path
        json_text = json.dumps(filled_json)

        # Create four jobs for each project.
        for state in [Job.ERROR, Job.QUEUED]:
            job = make_job(owner=owner, analysis=fastq_analysis, project=project, state=state, json_text=json_text)

    return


def init_users(sender, **kwargs):
    """
    Creates admin users if these are not present.
    """
    from engine.models import User
    logger.info("Setting up users")

    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            user = User(first_name=name, email=email, username=util.get_uuid(),
                        is_superuser=True, is_staff=True)
            user.set_password(settings.SECRET_KEY)
            user.save()
            logger.info(f"creating admin user: user.email={user.email}, user.id={user.id}")

    # Create a regular test user.
    test_buddy, new = User.objects.get_or_create(email="testbuddy@lvh.me")
    test_buddy.set_password("testbuddy@lvh.me")
    test_buddy.save()

    # create multiple 3 users and put them in groups.
    # assign a diffrent project to each group and allow one user to add others to groups.


    logger.info(f"creating user: {test_buddy.email}")


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


class EngineConfig(AppConfig):
    name = 'engine'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_site, sender=self)
        post_migrate.connect(init_users, sender=self)
        post_migrate.connect(init_proj, sender=self)
        #logger.debug("EngineConfig done")



