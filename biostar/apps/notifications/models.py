from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django import forms
from django.db.models import signals
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.utils.timezone import utc
from django.utils.translation import ugettext_lazy as _

# Inspired by django-messages at https://github.com/arneb/django-messages

class MessageManager(models.Manager):

    def inbox_for(self, user):
        "Returns all messages that were received by the given user"
        return self.filter(recipient=user)

    def outbox_for(self, user):
        "Returns all messages that were sent by the given user."
        return self.filter(sender=user)

# Create your models here.
class Message(models.Model):
    """
    A private message from user to user
    """

    body = models.TextField(_("Body"))
    sender = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='sent_messages', verbose_name=_("Sender"))
    subject = models.CharField(_("Subject"), max_length=120)

    recipient = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='received_messages', null=True, blank=True, verbose_name=_("Recipient"))
    parent_msg = models.ForeignKey('self', related_name='next_messages', null=True, blank=True, verbose_name=_("Parent message"))
    sent_at = models.DateTimeField(_("sent at"), null=False)
    read_at = models.DateTimeField(_("read at"), null=True, blank=True)

    objects = MessageManager()

    def __unicode__(self):
        return self.subject

    def save(self, **kwargs):

        if not self.id:
            self.sent_at = datetime.datetime.utcnow().replace(tzinfo=utc)
        super(Message, self).save(**kwargs)

from biostar.apps.posts.models import Post

class Subscription(models.Model):
    "Connects a post to a message"

    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name=_("User"))
    post = models.ForeignKey(Post, verbose_name=_("User"))
    creation_date = models.DateTimeField(_("Creation date"), null=False)

    @staticmethod
    def add_subscription(post, user):
        "Creates a subscription to a post"
        sub = Subscription(post=post, user=user)
        sub.creation_date = datetime.datetime.utcnow().replace(tzinfo=utc)
        sub.save()

def inbox_count_for(user):
    """
    returns the number of unread messages for the given user but does not
    mark them seen
    """
    return Message.objects.filter(recipient=user, read_at__isnull=True).count()

def init_subscription(sender, instance, created, *args, **kwargs):
    "Creates a subscription to a post"
    if created:
        Subscription.add_subscription(post=instance, user=instance.author)

signals.post_save.connect(init_subscription, sender=Post)