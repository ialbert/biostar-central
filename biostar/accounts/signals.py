import logging

from django.db.models.signals import post_save, pre_save
from django.dispatch import receiver

from biostar.accounts.models import Profile, User
from biostar.accounts import util, tasks

logger = logging.getLogger("engine")

@receiver(pre_save, sender=User)
def create_uuid(sender, instance, *args, **kwargs):
    # Generate a unique username if it does not exist.
    instance.username = instance.username or util.get_uuid(8)


@receiver(post_save, sender=User)
def create_profile(sender, instance, created, raw, using, **kwargs):

    if created:

        # Set the username a simpler username.
        username = f"{instance.pk}"

        # Fix uid clashes.
        if Profile.objects.filter(uid=username):
            username = util.get_uuid(8)
            logger.info(f"username clash for pk={instance.pk} new uid={username}")

        # Make sure staff users are also moderators.
        role = (
            Profile.MANAGER
            if (instance.is_staff or instance.is_superuser)
            else Profile.READER
        )

        # Create a user profile associate with the user.
        Profile.objects.create(
            user=instance, uid=username, name=instance.first_name, role=role
        )

        try:
            # Minimize failures on signup.
            # Send the welcome message.
            user_ids = [instance.pk]
            tasks.create_messages(user_ids=user_ids, template="messages/welcome.md")
        except Exception as exc:
            logger.error(f"{exc}")

    instance.profile.add_watched()
