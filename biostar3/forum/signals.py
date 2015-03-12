__author__ = 'ialbert'

import logging
from django.db.models.signals import post_save
from django.contrib.auth import get_user_model
from .mailer import EmailTemplate
from django.conf import settings
from django.utils.timezone import utc
from django.contrib.sites.models import Site

from datetime import datetime
from . import auth, models, mailer

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

        # Get or add a group subscription for the user.
        groupsub = auth.groupsub_get_or_create(user=instance.author, usergroup=instance.group)

        # Get or add a post subscription for the user
        postsub =  auth.postsub_get_or_create(user=instance.author, post=instance, groupsub=groupsub)

        # Add a message body for the new post.
        site = Site.objects.get_current()

        context = dict(post=instance, site=site,
                       slug=instance.group.domain)

        # This is the body of the message that gets created.
        em = mailer.EmailTemplate("post_created_message.html", data=context)
        body = models.MessageBody.objects.create(
            author=instance.author, subject=em.subj, content=em.text, html=em.html,
        )
        #em.send(to="localhost")

        sub_select = models.PostSub.objects.filter

        # ALl subscribers get local messages.
        def message_generator(mb):
            subs = sub_select(post=instance).select_related("user").all()
            for sub in subs:
                yield models.Message(user=sub.user, body=mb)

        # Bulk insert for all messages.
        models.Message.objects.bulk_create(message_generator(body),)

        # Some message preferences qualify for an email.
        email_targets = sub_select(post=instance, pref__in=settings.MESSAGE_EMAIL_PREFS).\
            select_related("user").all()

        # Collect all reply tokens

        #models.ReplyToken.objects.bulk_create(tokens, batch_size=100)


post_save.connect(user_update, sender=models.User)
post_save.connect(post_created, sender=models.Post)