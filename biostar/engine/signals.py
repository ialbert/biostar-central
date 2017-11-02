import logging, uuid
from .const import *
import os
from django.conf import settings
from django.core import management
import hjson as json
from django.core.files import File
from . import util
from biostar.tools import defaults

logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


def get_uuid(limit=None):
    return str(uuid.uuid4())[:limit]


def init_proj(sender, **kwargs):
    """
    Project 1 must exist. It will store existing analyses.
    We also create one job for each registered analysis.
    """
    from biostar.engine import models, auth

    # Get the first admin user.
    admin_user = models.User.objects.filter(is_superuser=True).first()

    #admin_group = models.Group.objects.get(name=settings.ADMIN_GROUP_NAME)

    # Demo projects enabled.


    if settings.DEMO_PROJECT_UID:
        uid = settings.DEMO_PROJECT_UID
        demo_project = models.Project.objects.filter(uid=uid).first()
        if not demo_project:
            auth.create_project(user=admin_user, uid=uid,
                                name=defaults.DEMO_PROJECT_NAME, summary=defaults.DEMO_PROJECT_SUMMARY,
                                text=defaults.DEMO_PROJECT_TEXT)

def init_users(sender, **kwargs):
    """
    Creates admin users and groups if needed.
    """
    from biostar.engine.models import User, Group

    # Get the admin group.
    group, created = Group.objects.get_or_create(name=settings.ADMIN_GROUP_NAME)

    logger.info("Setting up users")

    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            user = User(first_name=name, email=email,
                    is_superuser=True, is_staff=True)
            user.set_password(settings.SECRET_KEY)
            user.save()
            logger.info(f"Created admin user: {user.email}")
            group.user_set.add(user)

    if settings.DEBUG:
        testbuddy = 'testbuddy@lvh.me'
        user, flag = User.objects.get_or_create(email=testbuddy, username=testbuddy)
        user.set_password(testbuddy)
        user.save()

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
