from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver


@receiver(post_save, sender=User)
def finalize_post(sender, instance, created, **kwargs):
    if created:
        # Create subscriptions to the post using a task

        # Create notifications to the post using a task
        pass
