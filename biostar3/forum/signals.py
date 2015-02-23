__author__ = 'ialbert'

import logging
from django.contrib.auth.models import Permission
from django.db.models.signals import post_save
from django.contrib.auth import get_user_model
from .mailer import EmailTemplate
from . import models

User = models.User

logger = logging.getLogger("biostar")

# Moreator group created.
#

def user_update(sender, instance, created, **kwargs):

    if created:
        logger.info("created %s" % instance)
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