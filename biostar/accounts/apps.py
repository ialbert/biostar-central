import logging
import bleach

from django.apps import AppConfig
from django.conf import settings
from django.db.models.signals import post_migrate

logger = logging.getLogger('engine')


class AccountsConfig(AppConfig):
    name = 'biostar.accounts'

    def ready(self):
        from . import signals
        # Triggered after a migration command.
        post_migrate.connect(init_app, sender=self)


def init_app(sender, **kwargs):
    init_site()
    init_users()
    init_social()
    #init_tags()


def init_social():
    """Initialize social account providers."""

    from allauth.socialaccount.models import SocialApp
    from allauth.socialaccount.providers import registry
    from django.contrib.sites.models import Site

    # Populate the provider map based existing apps.
    providers = dict()
    for provider in registry.get_list():
        providers[provider.name] = provider

    # Create social apps as needed.
    for client in settings.SOCIAL_CLIENTS:

        name, client_id, client_secret = client

        # Check the app for existence.
        app = SocialApp.objects.filter(name=name)

        # Update the id and secrets to apply any changes that might have been made.
        if app.exists():
            SocialApp.objects.filter(name=name).update(client_id=client_id, secret=client_secret)
            continue

        # Create a new social app.
        logger.info(f"Creating social social app: {name}")

        # Identify the social app.
        provider = providers.get(name)

        if not provider:
            logger.error(f"Invalid provider name: {name}")
            continue

        # Create the provider here.
        site = Site.objects.filter(domain=settings.SITE_DOMAIN).first()

        app = SocialApp.objects.create(provider=provider.id, client_id=client_id, name=name,
                                       secret=client_secret)
        app.sites.add(site)
        app.save()


def init_users():
    """
    Creates admin users if needed.
    """
    from .models import User, Profile

    for name, email in settings.ADMINS:
        user = User.objects.filter(email=email).first()
        if not user:
            user = User.objects.create(first_name=name, email=email, is_superuser=True, is_staff=True)
            User.objects.filter(pk=user.pk).update(username=f'admin-{user.pk}')
            # Reload to update state that signals may change.
            user = User.objects.filter(pk=user.pk).first()
            user.set_password(settings.DEFAULT_ADMIN_PASSWORD)

            user.save()

            text = "Admin user created automatically on startup."
            Profile.objects.filter(user=user).update(name=name, text=text, html=text)
            logger.info(f"Creating admin user: {user.email}")
        else:

            User.objects.filter(pk=user.pk).update(is_superuser=True, is_staff=True)

            #Profile.objects.filter(user__pk=user.pk).update(state=Profile.TRUSTED)

            # You might want to reapply the default ADMIN password on migration.
            # This will log out the user from their current session.
            #user.set_password(settings.DEFAULT_ADMIN_PASSWORD)
            #user.save()

            #logger.info(f"Resetting password for admin user: {user.email}, {user.username}")
            logger.info(f"Admin user: {user.email} already exists")
            pass


def init_site():
    """
    Updates site domain and name.
    """
    from django.contrib.sites.models import Site

    # Print information on the database.
    db = settings.DATABASES['default']
    logger.info(f"db.name={db['NAME']}, db.engine={db['ENGINE']}")
    logger.info(f"email.backend={settings.EMAIL_BACKEND}, email.sender={settings.DEFAULT_FROM_EMAIL}")

    # Create the default site if necessary.
    Site.objects.get_or_create(id=settings.SITE_ID)

    # Update the default site domain and name.
    Site.objects.filter(id=settings.SITE_ID).update(domain=settings.SITE_DOMAIN, name=settings.SITE_NAME)

    # Get the current site
    site = Site.objects.get(id=settings.SITE_ID)
    logger.info("site.name={}, site.domain={}".format(site.name, site.domain))


def init_tags():
    """
    Initialize watched tags.

    Not possible to do in migrations:
    https://github.com/jazzband/django-taggit/issues/454
    """
    from .models import Profile

    profiles = Profile.objects.all()
    # TODO: add a progress to this.
    for pro in profiles:
        pro.add_watched()
        pro.save()

    logger.info("Added watched tags to profiles")
