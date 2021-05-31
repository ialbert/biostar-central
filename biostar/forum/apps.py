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
        post_migrate.connect(init_herald, sender=self)


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


def init_herald(sender, **kwargs):
    """
    Initialize the Biostar Herald Blog and Email Group.
    """
    from biostar.planet.models import Blog
    from django.shortcuts import reverse
    from biostar.emailer.models import EmailGroup

    title = "Biostar Herald"
    link = reverse('post_topic', kwargs=dict(topic='herald'))
    desc = "Share bioinformatics resources from across the web."

    uid = 'herald'

    hblog = Blog.objects.filter(link=link)
    group = EmailGroup.objects.filter(uid=uid)

    if hblog:
        hblog.update(desc=desc, title=title, remote=False)
    else:
        Blog.objects.create(title=title, desc=desc, remote=False, link=link)
        logger.info("created the biostar herald blog")

    if group:
        group.update(uid=uid, name=title, text=desc)
    else:
        EmailGroup.objects.create(uid=uid)
        logger.info("created the biostar herald email group")
