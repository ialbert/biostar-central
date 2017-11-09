from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings
import logging
from biostar.engine import util

logger = logging.getLogger('engine')

class AccountsConfig(AppConfig):
    name = 'biostar.accounts'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_site, sender=self)
        post_migrate.connect(init_users, sender=self)
        pass


def init_users(sender, **kwargs):
    """
    Creates admin users and groups if needed.
    """
    from biostar.engine.models import User, Group

    admin_group, created = Group.objects.get_or_create(name=settings.ADMIN_GROUP_NAME)

    logger.info("Setting up users")

    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            user = User(first_name=name, email=email,
                    is_superuser=True, is_staff=True)
            user.set_password(settings.SECRET_KEY)
            user.save()
            logger.info(f"Created admin user: {user.email}")
            admin_group.user_set.add(user)

    if settings.DEBUG:
        testbuddy = 'testbuddy@lvh.me'
        user, flag = User.objects.get_or_create(email=testbuddy, username=testbuddy)
        user.set_password(testbuddy)
        user.save()

    # Hardcoding a few users for now
    # TODO: move it to a command to add users
    users = [
        ('Aaron Maloy', 'aaron_maloy@fws.gov'),
        ('Meredith Bartron', 'meredith_bartron@fws.gov'),
        ('Doug Cavener', 'drc9@psu.edu'),
        ('Lan Wu Cavener', 'lxw34@psu.edu'),
        ]

    for name, email in users:
        user, flag = User.objects.get_or_create(email=email, username=util.get_uuid(8))
        if flag:
            user.set_password("testbuddy11")
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

