from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings
import logging

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
    Creates admin users if needed.
    """
    from .models import User

    logger.info("Setting up admin users")

    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            user = User.objects.create(first_name=name, email=email, is_superuser=True, is_staff=True)
            user.set_password(settings.DEFAULT_ADMIN_PASSWORD)
            user.save()
            logger.info(f"Created admin user: {user.email}")


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

