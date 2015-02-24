__author__ = 'ialbert'

import logging
from django.contrib.auth.models import Permission
from django.db.models.signals import post_save
from django.contrib.auth import get_user_model
from .mailer import EmailTemplate
from django.conf import settings
from django.utils.timezone import utc
from datetime import datetime

from . import models

User = models.User

logger = logging.getLogger("biostar")

def now():
    return datetime.utcnow().replace(tzinfo=utc)

def user_update(sender, instance, created, **kwargs):

    if created:
        logger.info("created %s" % instance)

        # Every user is a member of the default group.
        group = models.Group.objects.filter(name=settings.DEFAULT_GROUP_NAME).first()
        instance.groups.add(group)

        # Add a user profile on creation.
        right_now = now()
        profile = models.Profile.objects.create(
            user=instance, last_login=right_now, date_joined=right_now
        )

        # Send a welcome email to the user.
        data = dict(user=instance)
        em = EmailTemplate("user_creation.html", data=data)
        em.send(to=[instance.email])

    # Update moderator and admin group memberships on every save.
    mod_group, flag = models.get_or_create_group(name=models.MODERATOR_GROUP_NAME, user=instance)
    admin_group, flag = models.get_or_create_group(name=models.ADMIN_GROUP_NAME, user=instance)

    if instance.type == User.MODERATOR:
        instance.groups.add(mod_group)
    else:
        instance.groups.remove(mod_group)

    if instance.type == User.ADMIN:
        instance.groups.add(admin_group)
    else:
        instance.groups.remove(admin_group)

post_save.connect(user_update, sender=User)