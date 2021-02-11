import logging
from django.db.models.signals import post_migrate
from django.conf import settings
from django.apps import AppConfig
from django.template.defaultfilters import slugify

logger = logging.getLogger('engine')


class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        from . import signals
        # Triggered upon app initialization.
        post_migrate.connect(init_awards, sender=self)


def init_awards(sender, **kwargs):
    "Initializes the badges"
    from biostar.forum.models import Badge
    from biostar.forum.awards import ALL_AWARDS
    from biostar.accounts.models import Profile

    for obj in ALL_AWARDS:
        badge = Badge.objects.filter(name=obj.name)
        if badge:
            continue
        name = slugify(obj.name)
        badge = Badge.objects.create(name=obj.name, uid=name)

        # Badge descriptions may change.
        if badge.desc != obj.desc:
            badge.desc = obj.desc
            badge.icon = obj.icon
            badge.type = obj.type
            badge.save()

        logger.debug("initializing badge %s" % badge)
