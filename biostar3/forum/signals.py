__author__ = 'ialbert'

import logging
from django.db.models.signals import post_save
from django.contrib.auth import get_user_model
from .mailer import EmailTemplate
from django.conf import settings
from django.utils.timezone import utc
from django.contrib.sites.models import Site

from datetime import datetime
from . import auth, mailer
from models import *

logger = logging.getLogger("biostar")


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def user_update(sender, instance, created, **kwargs):
    if created:
        logger.info("created %s" % instance)

        # Every user will be a member of the default group.
        usergroup = UserGroup.objects.get(domain=settings.DEFAULT_GROUP_DOMAIN)

        # Create a subscription of the user to the default group.
        GroupSub.objects.create(user=instance, usergroup=usergroup)

        # Add a user profile on creation.
        right_now = now()
        profile = Profile.objects.create(
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

        # Default mode is Email for toplevel posts and local message otherwise.
        if instance.is_toplevel and groupsub.pref == settings.DEFAULT_MESSAGES:
            pref = settings.EMAIL_TRACKER
        else:
            pref = groupsub.pref

        # Get or add a post subscription for the user
        postsub = auth.postsub_get_or_create(user=instance.author, post=instance, pref=pref)

        # Add a message body for the new post.
        site = Site.objects.get_current()

        context = dict(post=instance, site=site,
                       slug=instance.group.domain)

        # This is the body of the message that gets created.
        em = mailer.EmailTemplate("post_created_message.html", data=context)
        body = MessageBody.objects.create(
            author=instance.author, subject=em.subj, content=em.text, html=em.html,
        )
        # em.send(to="localhost")

        sub_select = PostSub.objects.filter

        # ALl subscribers get local messages.
        def message_generator(mb):
            subs = sub_select(post=instance).select_related("user").all()
            for sub in subs:
                yield Message(user=sub.user, body=mb)

        # Bulk insert for all messages.
        Message.objects.bulk_create(message_generator(body), batch_size=100)

        # Some message preferences qualify for an email.
        def token_generator(obj):
            date = now()
            subs = sub_select(post=instance, pref__in=settings.MESSAGE_EMAIL_PREFS).select_related("user").all()
            for sub in subs:
                token = auth.make_uuid(size=8)
                em.send(to=[sub.user.email], token=token)
                yield ReplyToken(user=sub.user, post=obj, token=token, date=date)

        # Insert the reply tokens into the database.
        ReplyToken.objects.bulk_create(token_generator(instance), batch_size=100)


post_save.connect(user_update, sender=User)
post_save.connect(post_created, sender=Post)