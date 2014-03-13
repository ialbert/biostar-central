"""
This module is responsible for sending messages
to users based on various events
"""

from django.conf import settings
from biostar.apps.users.models import User, Profile
from biostar.apps.messages.models import Message, MessageBody
from biostar.apps.posts.models import Post, Tag
from django.core import mail
from biostar.const import *

# A basic message sending utility function.
def send_private_message(sender, recipient, content, sent_at=None):

    subject = "Private message from ?"
    sent_at = sent_at or now()
    body = MessageBody.objects.create(author=sender, subject=subject,
                                      text=content, sent_at=sent_at)

    msg = Message.objects.create(user=recipient, body=body, sent_at=body.sent_at)

    if recipient.profile.message_prefs == Profile.EMAIL_MESSAGE:
        try:
            # Bulk sending email messages.
            emails = [msg.email_tuple(recipient_list=recipient)]
            results = mail.send_mass_mail(emails)
        except Exception, exc:
            logger.error("email error %s" % exc)