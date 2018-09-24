import uuid
import mistune

from django.contrib.auth.models import User

from django.db.models.signals import post_save
from django.dispatch import receiver
from django.db import models
from biostar import settings


MAX_UID_LEN = 32
MAX_NAME_LEN = 80


def generate_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


class Manager(models.Manager):

    def get_all(self, **kwargs):
        "Return everything"
        return super().get_queryset().filter(**kwargs)


class Profile(models.Model):

    NEW, TRUSTED, SUSPENDED, BANNED = range(4)
    STATE_CHOICES = [(NEW, "New"), (TRUSTED, "Active"), (SUSPENDED, "Suspended"), (BANNED, "Banned")]
    state = models.IntegerField(default=NEW, choices=STATE_CHOICES, db_index=True)

    NORMAL, MODERATOR, MANAGER, BLOG = range(4)
    ROLE_CHOICES = [(NORMAL, "User"), (MODERATOR, "Moderator"), (MANAGER, "Manager"),
                    (BLOG, "Blog User")]

    NO_DIGEST, DAILY_DIGEST, WEEKLY_DIGEST, MONTHLY_DIGEST = range(4)
    DIGEST_CHOICES = [(NO_DIGEST, 'Never'), (DAILY_DIGEST, 'Daily'),
                      (WEEKLY_DIGEST, 'Weekly'), (MONTHLY_DIGEST, 'Monthly')]

    LOCAL_MESSAGE, EMAIL_MESSAGE, DIGEST_MESSAGES = range(3)
    MESSAGING_TYPE_CHOICES = [
                            (LOCAL_MESSAGE, "Local messages"),
                            (EMAIL_MESSAGE, "Email messages"),
                            (DIGEST_MESSAGES, "Digest email messages")
                            ]

    user = models.OneToOneField(User, on_delete=models.CASCADE)
    uid = models.CharField(max_length=MAX_UID_LEN, unique=True)
    name = models.CharField(max_length=MAX_NAME_LEN, default='', db_index=True)

    # Maximum amount of uploaded files a user is allowed to aggregate, in mega-bytes.
    max_upload_size = models.IntegerField(default=0)

    role = models.IntegerField(default=NORMAL, choices=ROLE_CHOICES)
    last_login = models.DateTimeField(null=True, db_index=True)

    # The number of new messages for the user.
    new_messages = models.IntegerField(default=0, db_index=True)

    # The last visit by the user.
    date_joined = models.DateTimeField(auto_now_add=True)

    # User provided location.
    location = models.CharField(default="", max_length=255, blank=True, db_index=True)

    # User provided website.
    website = models.URLField(default="", max_length=255, blank=True)

    # Google scholar ID
    scholar = models.CharField(default="", max_length=255, blank=True)

    score = models.IntegerField(default=0, db_index=True)

    # Twitter ID
    twitter = models.CharField(default="", max_length=255, blank=True)

    # This field is used to select content for the user.
    my_tags = models.CharField(default="", max_length=100, blank=True)

    # Description provided by the user html.
    text = models.TextField(default="", null=True, blank=True)

    html = models.TextField(null=True, blank=True)

    email_verified = models.BooleanField(default=False)

    # Notify when the a recipe has been complete.
    notify = models.BooleanField(default=False)
    message_prefs = models.IntegerField(choices=MESSAGING_TYPE_CHOICES,
                                        default=LOCAL_MESSAGE)

    # Subscription to daily and weekly digests.
    digest_prefs = models.IntegerField(choices=DIGEST_CHOICES, default=WEEKLY_DIGEST)

    # Opt-in to all messages from the site
    opt_in = models.BooleanField(default=False)

    objects = Manager()

    def __str__(self):
        return self.user.email

    def save(self, *args, **kwargs):
        self.uid = self.uid or generate_uuid(8)
        self.html = self.html or mistune.markdown(self.text)
        self.max_upload_size = self.max_upload_size or settings.MAX_UPLOAD_SIZE
        self.name = self.name or self.user.first_name or self.user.email.split("@")[0]
        super(Profile, self).save(*args, **kwargs)

    def can_moderate(self, source):
        "Check if the source user can moderate the target( self.user)"

        if source == self.user:
            return False

        return source.is_authenticated and (source.profile.is_manager or source.profile.is_moderator)

    @property
    def is_moderator(self):
        # Managers can moderate as well.
        return self.role == self.MODERATOR or self.role == self.MANAGER

    @property
    def is_manager(self):
        return self.role == self.MANAGER

    @property
    def is_suspended(self):
        return self.state == self.SUSPENDED


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):

    if created:
        # Make sure staff users are also moderators.
        role = Profile.MANAGER if instance.is_staff else Profile.NORMAL
        Profile.objects.create(user=instance, name=instance.first_name, role=role)

    instance.username = instance.username or f"user-{instance.pk}"