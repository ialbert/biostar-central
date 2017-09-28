from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate
import logging, uuid
import random


# Random block of latin text (lorem ipsum) to populate text feild
TEXT = """#TEST\nSed ut perspiciatis unde omnis iste natus error\n
sit voluptatem accusantium doloremque\nlaudantium, totam rem aperiam, eaque ipsa quae ab illo
"""

logger = logging.getLogger('engine')


def get_uuid(limit=None):
    return str(uuid.uuid4())[:limit]


def init_proj(sender, **kwargs):
    """
    Populate initial projects with N number data and analysis models
    """
    from engine.models import Project, Data, Analysis, Result
    from engine.models import User

    N = 2
    owner = User.objects.all().first()
    projects = [f"Project {x}" for x in range(1, 6)]
    data_analysis = [(f"Data {x}",f"Analysis {x}", f"Job {x}")
                     for x in range(0, 10)]

    for title in projects:

        test_set = random.sample(data_analysis, N)
        proj, flag = Project.objects.get_or_create(title=title, owner=owner, text=TEXT)
        proj.save()

        for data_title, analysis_title, job_title in test_set:

            datainput, flag = Data.objects.get_or_create(title=data_title,
                                                         owner=owner,
                                                         text=TEXT, project=proj)
            datainput.save()
            params, flag = Analysis.objects.get_or_create(title=analysis_title,
                                                          owner=owner,
                                                          text=TEXT, project=proj)
            params.save()

            # res, flag = Result.objects.get_or_create(title=job_title,
            #                                         owner=owner,
            #                                         text=TEXT, analysis=params)
            # res.save()

        logger.info(f'creating: {proj.title}')
    logger.info(f' with: {len(test_set)} data, analysis, and pipelines.')


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



