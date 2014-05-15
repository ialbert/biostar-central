'''
Inspired by django-messages at https://github.com/arneb/django-messages
'''

from __future__ import print_function, unicode_literals, absolute_import, division
import logging, datetime
from django.db import models
from django.conf import settings
from django.contrib import admin
from django.utils.timezone import utc
from django.utils.translation import ugettext_lazy as _
from django.core import mail

logger = logging.getLogger(__name__)

def now():
    return datetime.datetime.utcnow().replace(tzinfo=utc)

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
    MAX_SIZE = 120

    text = models.TextField(_("Text"))
    author = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='sent_messages', verbose_name=_("Sender"))
    subject = models.CharField(_("Subject"), max_length=MAX_SIZE)
    parent_msg = models.ForeignKey('self', related_name='next_messages', null=True, blank=True, verbose_name=_("Parent message"))
    sent_at = models.DateTimeField(_("sent at"), null=False)

    objects = MessageManager()

    def __unicode__(self):
        return self.subject

    def save(self, **kwargs):
        self.subject = self.subject[:self.MAX_SIZE]
        self.sent_at= self.sent_at or now()
        super(MessageBody, self).save(**kwargs)


# This contains the notification types.
from biostar.const import LOCAL_MESSAGE, MESSAGING_TYPE_CHOICES

# Connects user to message bodies
class Message(models.Model):
    "Connects recipents to sent messages"
    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='recipients', verbose_name=_("Recipient"))
    body = models.ForeignKey(MessageBody, related_name='messages', verbose_name=_("Message"))
    type = models.IntegerField(choices=MESSAGING_TYPE_CHOICES, default=LOCAL_MESSAGE, db_index=True)
    unread = models.BooleanField(default=True)
    sent_at = models.DateTimeField(db_index=True, null=True)

    def save(self, *args, **kwargs):
        self.sent_at = self.body.sent_at
        super(Message, self).save(**kwargs)

    def __unicode__(self):
        return u"Message %s, %s" % (self.user, self.body_id)

    @staticmethod
    def inbox_count_for(user):
        "Returns the number of unread messages for the given user but does not mark them seen"
        return MessageBody.objects.filter(recipient=user, unread=True).count()

    def email_tuple(self, recipient_list, from_email=None):
        "Returns an email tuple suitable to be mass emailed"
        from_email = from_email or settings.DEFAULT_FROM_EMAIL
        data = (self.body.subject, self.body.text, settings.DEFAULT_FROM_EMAIL, recipient_list)
        return data

# Admin interface to Message and MessageBody.
class MessageBodyAdmin(admin.ModelAdmin):
    search_fields = ('sender__name', 'sender__email', 'recipient__name', 'recipient__email', 'subject')
    list_select_related = ["sender", "post"]

# Admin interface to MessageBody
class MessageAdmin(admin.ModelAdmin):
    search_fields = ('recipient__name', 'recipient__email', 'recipient__name', 'recipient__email', 'subject')
    list_select_related = ["user", "post"]

#admin.site.register(Message, MessageAdmin)
admin.site.register(MessageBody, MessageBodyAdmin)


