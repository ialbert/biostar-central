import logging, uuid
import os
import hjson as json
from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate
from django.template.loader import get_template
from . import util
from engine.const import *

logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


def get_uuid(limit=None):
    return str(uuid.uuid4())[:limit]

def make_analysis_from_spec(path, user, project):
    from engine.models import Analysis
    json_data = open(path).read()
    json_obj = json.loads(json_data)
    title = json_obj["analysis_spec"]["title"]
    text = json_obj["analysis_spec"]["text"]
    template_path = json_obj["template"]["path"]
    makefile_template = get_template(template_path).template.source
    analysis = Analysis(json_data=json_data, owner=user, title=title, text=text,
             makefile_template=makefile_template, project=project)
    analysis.save()

    return analysis

def init_proj(sender, **kwargs):
    """
    Populate initial projects with N number data
    Creates one analysis model to allow for jobs to be run
    """
    from engine.models import Project, Data, Analysis, Job
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
                for state in (Job.RUNNING, Job.ERROR, Job.QUEUED):
                    title = analysis.title

                    job = Job(title=title, state=state,
                              project=project, analysis=analysis, owner=owner, makefile_template=analysis.makefile_template)
                    job.save()

    return


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



