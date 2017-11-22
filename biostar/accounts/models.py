from django.contrib.auth.models import Group
from django.contrib.auth.models import User

from django.db.models.signals import post_save
from django.dispatch import receiver
from django.db import models
from django.conf import settings
import uuid

MAX_UID_LEN = 32
MAX_NAME_LEN = 80

def generate_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

class Profile(models.Model):
    user = models.ForeignKey(User)
    uid = models.CharField(max_length=MAX_UID_LEN, unique=True)
    name = models.CharField(max_length=MAX_NAME_LEN, default='')

    def save(self, *args, **kwargs):
        self.uid = self.uid or generate_uuid(8)
        super(Profile, self).save(*args, **kwargs)

@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):

    if created:
        # Create a profile for user
        Profile.objects.create(user=instance, name=instance.first_name)


