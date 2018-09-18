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
        post_migrate.connect(init_social, sender=self)
        pass


def init_social(sender, **kwargs):
    """Initialize social account providers."""

    from allauth.socialaccount.models import SocialApp
    from allauth.socialaccount.providers.google.provider import GoogleProvider
    from django.contrib.sites.models import Site

    provider = GoogleProvider
    name = "google"

    client_id = settings.CLIENT_ID
    client_secret = settings.CLIENT_SECRET
    site = Site.objects.filter(domain=settings.SITE_DOMAIN)

    social_app = SocialApp.objects.filter(provider=provider.id, client_id=client_id,
                                          secret=client_secret, sites__in=site)
    if social_app.exists():
        return

    social_app = SocialApp.objects.create(provider=provider.id, client_id=client_id, name=name,
                                          secret=client_secret)
    social_app.sites.add(site.first())
    social_app.save()


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

