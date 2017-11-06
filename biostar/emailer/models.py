
from django.db import models
from biostar.accounts.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver
from biostar.engine.util import get_uuid



class EmailList(models.Model):

    ACTIVE, DELETED, INACTIVE, OPTOUT = 1,2,3,4

    STATE_CHOICES = [(ACTIVE, "Active"), (DELETED, "Deleted"), (INACTIVE, "Inactive"), (OPTOUT, "Opted out")]

    uid = models.CharField(max_length=16, blank=True, unique=True, default=get_uuid(16))
    state = models.IntegerField(default=ACTIVE, choices=STATE_CHOICES)



@receiver(post_save, sender=User)
def add_to_emaillist(sender, instance, created, **kwargs):

    if created:
        # Create a profile for user
        profile = Profile.objects.create(user=instance)