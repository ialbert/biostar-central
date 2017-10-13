import logging, uuid
import os
from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate
from engine.const import *

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
    from engine.models import Project, Data, make_analysis_from_spec, make_job, Job
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
            data = Data(title=data_title, owner=owner, text=data_desc, project=project, file=data_file, data_type=data_type)
            data.save()


    # Initialize the analyses.
    for test_spec in TEST_SPECS:
        analysis = make_analysis_from_spec(test_spec, user=owner, project=project)

        # Create four jobs for each project.
        for project in Project.objects.all():
            for state in (Job.RUNNING, Job.ERROR, Job.QUEUED, Job.FINISHED):
                # user, analysis, and project are the only necessary things
                job = make_job(owner=owner, analysis=analysis, project=project, state=state)
                job.save()
                #logger.info(f'Created job: {job.title} in project : {project.title} with state : {job.get_state_display()}')

    return


def init_users(sender, **kwargs):
    """
    Creates admin users if these are not present.
    """
    from engine.models import User
    logger.info("Setting up users")

    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            user = User(first_name=name, email=email, username=get_uuid(16),
                        is_superuser=True, is_staff=True)
            user.set_password(settings.SECRET_KEY)
            user.save()
            logger.info(f"creating admin user: user.email={user.email}, user.id={user.id}")

    # Create a regular test user.
    #testbuddy @ lvh.me
    test_buddy, new = User.objects.get_or_create(email="testbuddy@lvh.me")
    test_buddy.set_password("testbuddy@lvh.me")
    test_buddy.save()

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
        logger.debug("EngineConfig done")



