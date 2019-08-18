import random

from django.apps import AppConfig
from django.conf import settings

from django.db.models.signals import post_migrate


class EmailerConfig(AppConfig):
    name = 'biostar.emailer'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init, sender=self)


def init(sender, **kwargs):
    """
    Initialize with groups.
    """
    from .models import EmailAddress, EmailGroup, Subscription, get_uuid
    # Generate random groups and subscriptions

    for num in range(1, 5):
        name = f"Group-{num}"
        group, flag = EmailGroup.objects.get_or_create(name=name, text=name, html=name)

    for num in range(1, 20):
        uid = get_uuid(8)
        name, email = f"Name-{num}", f"email-{uid}j@nomail.for.me"
        address, flag = EmailAddress.objects.get_or_create(name=name, email=email)

    groups = list(EmailGroup.objects.all()[:20])
    addresses = list(EmailAddress.objects.all()[:20])

    for num in range(1, 20):
        group = random.choice(groups)
        address = random.choice(addresses)
        sub, flag = Subscription.objects.get_or_create(address=address, group=group)
