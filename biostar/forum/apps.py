import logging
from django.db.models.signals import post_migrate
from django.conf import settings
from django.apps import AppConfig

logger = logging.getLogger('engine')


class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        from . import signals
        # Triggered upon app initialization.
        post_migrate.connect(init_awards, sender=self)
        post_migrate.connect(init_digest, sender=self)


def init_digest(sender, **kwargs):
    from biostar.accounts.models import Profile

    # Ensure digest emails are not sent when debugging
    if settings.DEBUG:
        profiles = Profile.objects.all()
        Profile.objects.filter(id__in=profiles).update(digest_prefs=Profile.NO_DIGEST)


def init_awards(sender, **kwargs):
    "Initializes the badges"
    from biostar.forum.models import Badge
    from biostar.forum.awards import ALL_AWARDS
    from biostar.accounts.models import Profile

    for obj in ALL_AWARDS:
        badge = Badge.objects.filter(name=obj.name)
        if badge:
            continue
        badge = Badge.objects.create(name=obj.name)

        # Badge descriptions may change.
        if badge.desc != obj.desc:
            badge.desc = obj.desc
            badge.icon = obj.icon
            badge.type = obj.type
            badge.save()

        logger.debug("initializing badge %s" % badge)
