from django.db.models.signals import post_save, pre_save
from django.dispatch import receiver

from biostar.accounts.models import Profile, User
from biostar.accounts import util, tasks


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, raw, using, **kwargs):
    if created:
        # Set the username to a simpler form.
        username = f"user-{instance.pk}"
        User.objects.filter(pk=instance.pk).update(username=username)

        # Make sure staff users are also moderators.
        role = Profile.MANAGER if instance.is_staff else Profile.READER
        Profile.objects.using(using).create(user=instance, uid=username, name=instance.first_name, role=role)
        tasks.create_messages(rec_list=[instance], template="messages/welcome.html")


@receiver(pre_save, sender=User)
def create_uuid(sender, instance, *args, **kwargs):
    instance.username = instance.username or util.get_uuid(8)


