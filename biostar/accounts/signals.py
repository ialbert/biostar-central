from django.db.models.signals import post_save, pre_save
from django.dispatch import receiver

from biostar.accounts.models import Profile, User
from biostar.accounts import util, tasks
from django.template.defaultfilters import slugify


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, raw, using, **kwargs):
    if created:
        # Set the username to a simpler form.
        username = f"{instance.first_name}-{instance.pk}" if instance.first_name else f'user-{instance.pk}'
        username = slugify(username)
        if User.objects.filter(username=username).exclude(id=instance.pk).exists():
            username = util.get_uuid(6)

        User.objects.filter(pk=instance.pk).update(username=username)

        # Make sure staff users are also moderators.
        role = Profile.MANAGER if instance.is_staff else Profile.READER
        Profile.objects.using(using).create(user=instance, uid=username, name=instance.first_name, role=role)
        user_ids = [instance.pk]
        tasks.create_messages.spool(user_ids=user_ids, template="messages/welcome.md")

    # Recompute watched tags
    instance.profile.add_watched()


@receiver(pre_save, sender=User)
def create_uuid(sender, instance, *args, **kwargs):

    instance.username = instance.username or util.get_uuid(8)