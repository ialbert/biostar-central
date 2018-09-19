import logging

from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate

logger = logging.getLogger('engine')


class AccountsConfig(AppConfig):
    name = 'biostar.accounts'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_site, sender=self)
        post_migrate.connect(init_users, sender=self)
        post_migrate.connect(init_social, sender=self)
        pass


def init_social(sender, **kwargs):
    """Initialize social account providers."""

    from allauth.socialaccount.models import SocialApp
    from allauth.socialaccount.providers import registry
    from django.contrib.sites.models import Site

    # Populate the provider map based existing apps.
    providers = dict()
    for provider in registry.get_list():
        providers[provider.name.lower()] = provider

    # Create social apps as needed.
    for client in settings.SOCIAL_CLIENTS:

        name, client_id, client_secret = client

        app = SocialApp.objects.filter(client_id=client_id)

        # If app exists we are done.
        if app.exists():
            continue

        # Create this social app.
        logger.info(f"Creating social social app: {name}")

        # Identify the social app
        provider = providers.get(name.lower())

        if not provider:
            logger.error(f"Invalid provider name: {name}")
            continue

        # Create the provider here.
        site = Site.objects.filter(domain=settings.SITE_DOMAIN).first()

        app = SocialApp.objects.create(provider=provider.id, client_id=client_id, name=name,
                                       secret=client_secret)
        app.sites.add(site)
        app.save()


def init_users(sender, **kwargs):
    """
    Creates admin users if needed.
    """
    from .models import User, Profile

    logger.info("Setting up admin users")

    for name, email in settings.ADMINS:
        if not User.objects.filter(email=email):
            user = User.objects.create(first_name=name, email=email, is_superuser=True, is_staff=True)
            user.set_password(settings.DEFAULT_ADMIN_PASSWORD)
            user.save()
            text = "Local user started with the website"
            Profile.objects.filter(user__pk=user.pk).update(location="State College",
                                                            text=text,
                                                            html=text)
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
