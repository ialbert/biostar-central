__author__ = 'ialbert'

import logging
from django.db.models.signals import post_save
from django.contrib.auth import get_user_model
from .mailer import EmailTemplate
from django.conf import settings
from django.utils.timezone import utc
from datetime import datetime
from . import auth, models, mailer


User = models.User

logger = logging.getLogger("biostar")


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def user_update(sender, instance, created, **kwargs):
    if created:
        logger.info("created %s" % instance)

        # Every user will be a member of the default group.
        usergroup = models.UserGroup.objects.get(domain=settings.DEFAULT_GROUP_DOMAIN)

        # Create a subscription of the user to the default group.
        models.GroupSub.objects.create(user=instance, usergroup=usergroup)

        # Add a user profile on creation.
        right_now = now()
        profile = models.Profile.objects.create(
            user=instance, last_login=right_now, date_joined=right_now
        )

        if settings.SEND_WELCOME_EMAIL:
            # Send a welcome email to the user.
            data = dict(user=instance)
            em = EmailTemplate("user_welcome_email.html", data=data)
            em.send(to=[instance.email])


def post_created(sender, instance, created, **kwargs):
    # This is where messages are sent
    if created:
        logger.info("created %s" % instance)
        if instance.is_toplevel:
            auth.add_groupsub(user=instance.author, usergroup=instance.group)

        # Add a message body for the new post.
        context = dict(post=instance, slug=instance.group.domain)

        em = mailer.EmailTemplate("post_created_message.html", data=context)
        body = models.MessageBody.objects.create(
            author=instance.author,
            subject=em.subj,
            content=em.text,
            html=em.html,
        )

post_save.connect(user_update, sender=User)
post_save.connect(post_created, sender=models.Post)