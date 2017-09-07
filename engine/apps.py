from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate, pre_migrate
import logging, uuid

logger = logging.getLogger('engine')


def get_uuid(limit=None):
    return str(uuid.uuid4())[:limit]


def init_proj(sender, **kwargs):
    """
    Populate initial projects
    """
    from engine.models import User
    from engine.models import Project
    titles = ['Project 1', 'Project 2']

    owner = User.objects.all().first()
    for title in titles:
        proj, flag = Project.objects.get_or_create(title=title, owner=owner)
        logger.info(f'creating: {proj.title}')

def init_users(sender, **kwargs):
    """
    Creates admin users if these are not present.
    """
    from engine.models import User

    # Create the super users according to the settings.
    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            logger.info(f"creating admin user.name={user.name}, user.email={user.email}")
            user = User(email=email, username=get_uuid(16), is_superuser=True, is_staff=True)
            user.set_password(settings.SECRET_KEY)
            user.save()

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



