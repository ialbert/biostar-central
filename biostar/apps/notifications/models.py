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

# A message body is information sent to users.
class MessageBody(models.Model):
    """
    A private message from user to user
    """

    body = models.TextField(_("Body"))
    sender = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='sent_messages', verbose_name=_("Sender"))
    subject = models.CharField(_("Subject"), max_length=120)

    parent_msg = models.ForeignKey('self', related_name='next_messages', null=True, blank=True, verbose_name=_("Parent message"))
    sent_at = models.DateTimeField(_("sent at"), null=False)

    objects = MessageManager()

    def __unicode__(self):
        return self.subject

    def save(self, **kwargs):

        if not self.id:
            self.sent_at = datetime.datetime.utcnow().replace(tzinfo=utc)
        super(Message, self).save(**kwargs)


# Connects user to message bodies
class Message(models.Model):
    "Connects recipents to sent messages"
    recipient = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='recipients', verbose_name=_("Recipient"))
    message = models.ForeignKey(MessageBody, related_name='messages', verbose_name=_("Message"))
    read_at = models.DateTimeField(_("read at"), null=True, blank=True, db_index=True)

    @staticmethod
    def inbox_count_for(user):
        "Returns the number of unread messages for the given user but does not mark them seen"
        return MessageBody.objects.filter(recipient=user, read_at__isnull=True).count()

# Admin interface to Message and MessageBody.
class MessageBodyAdmin(admin.ModelAdmin):
    search_fields = ('sender__name', 'sender__email', 'recipient__name', 'recipient__email', 'subject')
    list_select_related = ["sender", "post"]

# Admin interface to MessageBody
class MessageAdmin(admin.ModelAdmin):
    search_fields = ('recipient__name', 'recipient__email', 'recipient__name', 'recipient__email', 'subject')
    list_select_related = ["user", "post"]

admin.site.register(Message, MessageAdmin)
admin.site.register(MessageBody, MessageBodyAdmin)

# Set up signals on a new post.
from biostar.apps.posts.models import Post
from . import actions

# Add the notification related signals
signals.post_save.connect(actions.new_post_created, sender=Post)



