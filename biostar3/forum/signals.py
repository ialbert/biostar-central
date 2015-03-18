__author__ = 'ialbert'

import logging
from django.db.models.signals import post_save
from django.contrib.auth import get_user_model
from .mailer import EmailTemplate
from django.conf import settings
from django.utils.timezone import utc
from django.contrib.sites.models import Site

from datetime import datetime
from . import auth, mailer, tasks
from models import *
from allauth.account.signals import user_logged_in

logger = logging.getLogger("biostar")

def user_login(sender, request, user, **kwargs):
    # Actions performed on user login
    tasks.add_user_location.delay(request, user)


def user_create(sender, instance, created, **kwargs):
    if created:
        logger.info("created %s" % instance)

        # Every user will be a member of the default group.
        usergroup = UserGroup.objects.get(domain=settings.DEFAULT_GROUP_DOMAIN)

        # Create a subscription of the user to the default group.
        GroupSub.objects.create(user=instance, usergroup=usergroup)
        now = right_now()
        # Add a user profile on creation.
        profile = Profile.objects.create(
            user=instance, last_login=now, date_joined=now,
        )

        if settings.SEND_WELCOME_EMAIL:
            # Send a welcome email to the user.
            data = dict(user=instance)
            em = EmailTemplate("user_welcome_email.html", data=data)
            em.send_email(to=[instance.email])


def post_created(sender, instance, created, **kwargs):
    # This is where messages are sent
    if created:
        logger.info("created %s" % instance)

        # Subscriptions will apply relative to the root.
        # Get or add the group subscription for the user.
        groupsub = auth.groupsub_get_or_create(user=instance.author, usergroup=instance.root.usergroup)

        # Get or add the post subscription for the user.
        postsub = auth.postsub_get_or_create(user=instance.author, post=instance.root, sub_type=groupsub.type)

        # Route the message creation via celery if necessary.
        if settings.CELERY_ENABLED:
            tasks.create_messages.delay(instance)
        else:
            tasks.create_messages(instance)

post_save.connect(user_create, sender=User)
post_save.connect(post_created, sender=Post)
user_logged_in.connect(user_login)