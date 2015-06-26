__author__ = 'ialbert'

import logging
from django.db.models.signals import post_save
from django.contrib.auth import get_user_model
from .mailer import EmailTemplate
from django.conf import settings
from django.utils.timezone import utc
from django.contrib.sites.models import Site
from django.db.models import Q
from datetime import datetime
from . import auth, mailer, tasks
from .models import *
from allauth.account.signals import user_logged_in
from django.dispatch import receiver

logger = logging.getLogger("biostar")


@receiver(user_logged_in)
def user_login(sender, request, user, **kwargs):
    # Actions performed on user login
    ip = auth.remote_ip(request)
    func = tasks.add_user_location
    func.delay(ip, user) if settings.CELERY_ENABLED else func(ip, user)


def user_create(sender, instance, created, **kwargs):
    if created:

        logger.info("%s" % instance)

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


def post_saved(sender, instance, created, **kwargs):

    if created and instance.is_toplevel:
        # Top level posts need to fill the root and parent ids.
        # Self referential ForeignKeys will not be set otherwise.
        Post.objects.filter(pk=instance.pk).update(root_id=instance.pk, parent_id=instance.pk)

    # Update the subscriptions upon post save.
    if settings.CELERY_ENABLED:
        tasks.update_post_subscriptions.delay(instance)
        tasks.send_notifications.delay(instance)
    else:
        tasks.update_post_subscriptions(instance)
        tasks.send_notifications(instance)

    # This is where messages are sent
    if created:
        logger.info("%s" % instance)

        if not instance.uuid:
            # If the unique id not set then set it to the primary key.
            Post.objects.filter(pk=instance.pk).update(uuid=instance.pk)

    # Update the reply count on the post.
    instance.update_reply_count()

def register():
    post_save.connect(user_create, sender=User, dispatch_uid="user_create")
    post_save.connect(post_saved, sender=Post, dispatch_uid="post_saved")
