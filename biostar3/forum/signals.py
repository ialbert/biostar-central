__author__ = 'ialbert'

import logging
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
        group = models.UserGroup.objects.filter(name=settings.DEFAULT_GROUP_NAME).first()
        instance.usergroups.add(group)

        # Add a user profile on creation.
        right_now = now()
        profile = models.Profile.objects.create(
            user=instance, last_login=right_now, date_joined=right_now
        )

        # Send a welcome email to the user.
        data = dict(user=instance)
        em = EmailTemplate("user_creation.html", data=data)
        em.send(to=[instance.email])

post_save.connect(user_update, sender=User)