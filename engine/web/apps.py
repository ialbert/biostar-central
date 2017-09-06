from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate
import logging, uuid

logger = logging.getLogger('engine')


'''
 if not user:
        # This user needs to be created/
        created = True
        user = User.objects.create(email=email, username=get_uuid())
        user.set_password(get_uuid())
        user.save()
        user.profile.name = name
'''

def get_uuid(limit=None):
    return str(uuid.uuid4())[:limit]


def init_users(sender, **kwargs):
    """
    Creates admin users if these are not present.
    """
    from engine.web.models import User

    # Create the super users according to the settings.
    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            logger.info("creating admin user.name={}, user.email={}".format(name, email))
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
    from engine.web.models import *

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
