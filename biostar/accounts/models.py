import uuid
from datetime import datetime, timedelta

import mistune
from django.conf import settings
from django.shortcuts import reverse
from django.contrib.auth.models import User
from django.db import models
from django.utils.timezone import utc
from biostar.accounts import util


def fixcase(name):
    return name.upper() if len(name) == 1 else name.lower()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


MAX_UID_LEN = 255
MAX_NAME_LEN = 255
MAX_TEXT_LEN = 10000
MAX_FIELD_LEN = 1024


class Profile(models.Model):
    NEW, TRUSTED, SUSPENDED, BANNED = range(4)
    STATE_CHOICES = [(NEW, "New"), (TRUSTED, "Active"), (SUSPENDED, "Suspended"), (BANNED, "Banned")]
    state = models.IntegerField(default=NEW, choices=STATE_CHOICES, db_index=True)

    READER, MODERATOR, MANAGER, BLOG = range(4)
    ROLE_CHOICES = [(READER, "Reader"), (MODERATOR, "Moderator"), (MANAGER, "Admin"), (BLOG, "Blog User")]

    NO_DIGEST, DAILY_DIGEST, WEEKLY_DIGEST, MONTHLY_DIGEST,ALL_MESSAGES  = range(5)

    DIGEST_CHOICES = [(NO_DIGEST, 'Never'), (DAILY_DIGEST, 'Daily'),
                      (WEEKLY_DIGEST, 'Weekly'), (MONTHLY_DIGEST, 'Monthly'),
                      (ALL_MESSAGES, "Email for every new thread (mailing list mode)")
                      ]

    LOCAL_MESSAGE, EMAIL_MESSAGE, NO_MESSAGES, DEFAULT_MESSAGES = range(4)
    MESSAGING_TYPE_CHOICES = [
        (DEFAULT_MESSAGES, "Default"),
        (EMAIL_MESSAGE, "Email"),
        (LOCAL_MESSAGE, "Local Messages"),
        (NO_MESSAGES, "No messages"),

    ]

    # Connection to the user.
    user = models.OneToOneField(User, on_delete=models.CASCADE)

    # User unique id (handle)
    uid = models.CharField(max_length=MAX_UID_LEN, unique=True)

    # User diplay name.
    name = models.CharField(max_length=MAX_NAME_LEN, default='', db_index=True)

    # Maximum amount of uploaded files a user is allowed to aggregate, in mega-bytes.
    max_upload_size = models.IntegerField(default=0)

    # The role of the user.
    role = models.IntegerField(default=READER, choices=ROLE_CHOICES)

    # The date the user last logged in.
    last_login = models.DateTimeField(null=True, max_length=255, db_index=True)

    # The number of new messages for the user.
    new_messages = models.IntegerField(default=0, db_index=True)

    # The last visit by the user.
    date_joined = models.DateTimeField(auto_now_add=True, max_length=255)

    # User provided location.
    location = models.CharField(default="", max_length=255, blank=True, db_index=True)

    # User provided website.
    website = models.URLField(default="", max_length=255, blank=True)

    # Google scholar ID
    scholar = models.CharField(default="", max_length=255, blank=True)

    # User reputation score.
    score = models.IntegerField(default=0, db_index=True)

    # Twitter ID
    twitter = models.CharField(default="", max_length=255, blank=True)

    # This field is used to select content for the user.
    my_tags = models.CharField(default="", max_length=255, blank=True)

    # Description provided by the user html.
    text = models.TextField(default="No profile information", null=True, max_length=MAX_TEXT_LEN, blank=True)

    # The html version of the user information.
    html = models.TextField(null=True, max_length=MAX_TEXT_LEN, blank=True)

    # The state of the user email verfication.
    email_verified = models.BooleanField(default=False)

    # Automatic notification
    notify = models.BooleanField(default=False)

    # Default subscriptions inherit from this
    message_prefs = models.IntegerField(choices=MESSAGING_TYPE_CHOICES, default=DEFAULT_MESSAGES)

    # Subscription to daily and weekly digests.
    digest_prefs = models.IntegerField(choices=DIGEST_CHOICES, default=WEEKLY_DIGEST)

    # Opt-in to all messages from the site
    opt_in = models.BooleanField(default=False)

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        self.uid = self.uid or util.get_uuid(8)
        self.html = self.html or mistune.markdown(self.text)
        self.max_upload_size = self.max_upload_size or settings.MAX_UPLOAD_SIZE
        self.name = self.name or self.user.first_name or self.user.email.split("@")[0]
        self.date_joined = self.date_joined or now()
        self.last_login = self.last_login or now() - timedelta(days=1)
        super(Profile, self).save(*args, **kwargs)

    @property
    def is_moderator(self):
        # Managers can moderate as well.
        return self.role == self.MODERATOR or self.role == self.MANAGER or self.user.is_staff or self.user.is_superuser

    @property
    def trusted(self):
        return self.user.is_staff or self.state == self.TRUSTED

    @property
    def is_manager(self):
        return self.role == self.MANAGER

    @property
    def get_absolute_url(self):

        return reverse('user_profile', kwargs=dict(uid=self.uid))

    @property
    def is_suspended(self):
        return self.state == self.SUSPENDED

# Connects user to message bodies
class MessageBody(models.Model):
    """
    A message that may be shared across all users.
    """
    body = models.TextField(max_length=MAX_TEXT_LEN)
    html = models.TextField(default='', max_length=MAX_TEXT_LEN * 10)

    def save(self, *args, **kwargs):
        self.html = self.html or mistune.markdown(self.body)
        super(MessageBody, self).save(**kwargs)

# Connects user to message bodies
class Message(models.Model):
    """
    Connects recipients to sent messages
    """

    SPAM, VALID, UNKNOWN = range(3)
    SPAM_CHOICES = [(SPAM, "Spam"), (VALID, "Not spam"), (UNKNOWN, "Unknown")]
    spam = models.IntegerField(choices=SPAM_CHOICES, default=UNKNOWN)

    uid = models.CharField(max_length=32, unique=True)
    sender = models.ForeignKey(settings.AUTH_USER_MODEL, related_name="author", on_delete=models.CASCADE)
    recipient = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE)

    subject = models.CharField(max_length=120)
    body = models.ForeignKey(MessageBody, on_delete=models.CASCADE)

    unread = models.BooleanField(default=True)
    sent_date = models.DateTimeField(db_index=True, null=True)

    def save(self, *args, **kwargs):
        self.uid = self.uid or util.get_uuid(10)
        self.sent_date = self.sent_date or util.now()
        super(Message, self).save(**kwargs)

    def __str__(self):
        return f"Message {self.sender}, {self.recipient}"
