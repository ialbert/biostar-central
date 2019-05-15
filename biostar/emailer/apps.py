import random
from django.db.models.signals import post_migrate
from django.apps import AppConfig


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

    for i in range(1, 5):
        name = f"Group-{i}"
        flag, group = EmailGroup.objects.get_or_create(name=name, text=name, html=name)

    groups = list(EmailGroup.objects.all()[:20])

    for j in range(1, 10):
        uid = get_uuid(4)
        group = random.choice(groups)
        name, email = f"Name-{j}", f"email-{uid}j@nomail.for.me"
        flag, email = EmailAddress.objects.get_or_create(name=name, email=email, group=group)


