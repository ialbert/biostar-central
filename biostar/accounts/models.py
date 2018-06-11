from django.contrib.auth.models import User

from django.db.models.signals import post_save
from django.dispatch import receiver
from django.db import models
from biostar import settings
import uuid

MAX_UID_LEN = 32
MAX_NAME_LEN = 80

def generate_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


class Profile(models.Model):

    user = models.OneToOneField(User, on_delete=models.CASCADE)
    uid = models.CharField(max_length=MAX_UID_LEN, unique=True)
    name = models.CharField(max_length=MAX_NAME_LEN, default='')

    # Maximum amount of uploaded files a user is allowed to aggregate, in mega-bytes.
    max_upload_size = models.IntegerField(default=0)

    NEW, TRUSTED, SUSPENDED, BANNED = range(1,5)
    STATE_CHOICES = [(NEW, "New"), (TRUSTED, "Active"), (SUSPENDED, "Suspended"), (BANNED, "Banned")]
    state = models.IntegerField(default=NEW, choices=STATE_CHOICES)

    NORMAL, MODERATOR, MANAGER = 1,2,3
    ROLE_CHOICES = [(NORMAL, "Normal User"),(MODERATOR, "Moderator"), (MANAGER, "Manager")]
    role = models.IntegerField(default=NORMAL, choices=ROLE_CHOICES)
    last_login = models.DateTimeField(auto_now_add=True)

    notify = models.BooleanField(default=False)

    def __str__(self):
        return self.user.email

    def save(self, *args, **kwargs):
        self.uid = self.uid or generate_uuid(8)
        self.max_upload_size = self.max_upload_size or settings.MAX_UPLOAD_SIZE
        self.name = self.user.first_name
        super(Profile, self).save(*args, **kwargs)

    @property
    def is_moderator(self):
        return self.role == self.MODERATOR

    @property
    def is_manager(self):
        return self.role == self.MANAGER



@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):

    if created:
        # Make sure staff users are also moderators.
        role = Profile.MANAGER if instance.is_staff else Profile.NORMAL
        Profile.objects.create(user=instance, name=instance.first_name, role=role)

    instance.username = instance.username or f"user-{instance.pk}"