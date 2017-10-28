import logging, uuid
from engine.const import *
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
    from engine.models import User, Group, Project

    # Get the first admin user.
    admin_user = User.objects.filter(is_superuser=True).first()
    admin_group = Group.objects.get(name=settings.ADMIN_GROUP_NAME)

    # The first project is admin project.
    admin_project = Project.objects.filter(id=1).first()

    # The second project is a demo project.
    demo_project = Project.objects.filter(id=2).first()

    # Set it up first time around.
    if not admin_project:
        admin_project = Project.objects.create(
            name=defaults.ADMIN_PROJECT_NAME,
            summary=defaults.ADMIN_PROJECT_SUMMARY,
            text=defaults.ADMIN_PROJECT_TEXT,
            owner=admin_user,
            type=Project.ADMIN
        )

    if not demo_project:
        demo_project = Project.objects.create(
            name=defaults.DEMO_PROJECT_NAME,
            summary=defaults.DEMO_PROJECT_SUMMARY,
            text=defaults.DEMO_PROJECT_TEXT,
            owner=admin_user,
            type=Project.USER
        )


def init_users(sender, **kwargs):
    """
    Creates admin users and groups if needed.
    """
    from engine.models import User, Group

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

    # TODO: remove later
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
