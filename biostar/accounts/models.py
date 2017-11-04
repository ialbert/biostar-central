from django.contrib.auth.models import Group
from django.contrib.auth.models import User

from django.db.models.signals import post_save
from django.dispatch import receiver
from django.db import models
from django.conf import settings
import uuid

def generate_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

class Profile(models.Model):
    user = models.ForeignKey(User)
    uid = models.CharField(max_length=32, unique=True)

    def save(self, *args, **kwargs):
        self.uid = self.uid or generate_uuid(8)
        super(Profile, self).save(*args, **kwargs)

@receiver(post_save, sender=User)
def create_profile(sender, instance, created, **kwargs):

    if created:
        # Create a profile for user
        profile = Profile.objects.create(user=instance)

        # Every user belongs to their own USER group.
        user_group, created = Group.objects.get_or_create(name=profile.uid)
        user_group.user_set.add(instance)

post_save.connect(create_profile, sender=User)
