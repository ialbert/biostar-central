from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate
import logging, uuid
from .settings import BASE_DIR
import os


# This is a temporary data structure.
TEST_PROJECTS = [
    # title, info pairs
    ("Sequencing run 1", "Lamar sequencing center"),
    ("Sequencing run 2", "Lamar sequencing center"),
    ("Sequencing run 3", "Lamar sequencing center"),
]


TEST_DATA = [
    ("Compressed data directory", "This directory contains all datasets for the run"),
    ("Sample sheet", "This file contains a sample sheet describing the data in the directory"),
]



logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


JSON_SPECFILE =join(BASE_DIR, '..', 'pipeline',
                'templates','qc', 'qc_spec.hjson' )


def get_uuid(limit=None):
    return str(uuid.uuid4())[:limit]


def init_proj(sender, **kwargs):
    """
    Populate initial projects with N number data
    Creates one analysis model to allow for jobs to be run
    """
    from engine.models import Project, Data, Analysis, Job
    from engine.models import User

    owner = User.objects.all().first()

    # Make a project
    for title, description in TEST_PROJECTS:

        project, flag = Project.objects.get_or_create(title=title, owner=owner, text=description)

        # add some data to the project
        for data_title, data_desc in TEST_DATA:
            data, flag = Data.objects.get_or_create(title=data_title,
                                                         owner=owner,
                                                         text=data_desc, project=project)
            data.save()

        logger.info(f'creating or getting: {project.title}')

    analysis, flag = Analysis.objects.get_or_create(title="Analysis 1",
                                                    owner=owner,
                                                    text="analysis description",
                                                    spec_origin=JSON_SPECFILE)
    #print(analysis.spec_source)
    analysis.save()
    #print(analysis.spec_source)

    # Pick most recent project to make a job out of
    jproject = Project.objects.order_by("-id").first()

    # Make a job in each state( 3 jobs to one project and analysis)
    states = {"Queued":1, "Running":2, "Finished":3, "Error":4}

    for state in states:

        job, flag = Job.objects.get_or_create(title="Result 1",
                                        text="job description",
                                        project=jproject,
                                        analysis=analysis,
                                        owner=owner,
                                        state=states[state])
        job.save()
        logger.info(f' job={job.id} made in {state} state')


def init_users(sender, **kwargs):
    """
    Creates admin users if these are not present.
    """
    from engine.models import User
    logger.info("Setting up users")

    for name, email in settings.ADMINS:
        
        if not User.objects.filter(email=email):
            user = User(first_name="foo", email=email, username=get_uuid(16),
                        is_superuser=True, is_staff=True)
            user.set_password(settings.SECRET_KEY)
            user.save()
            logger.info(f"creating admin: user.username = {user.username}, user.email={user.email}")

    # Create a regular test user.
    test_user = User.objects.get_or_create(email="foo@bar.com", password="foobar221")


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
        logger.debug("EngineConfig done")



