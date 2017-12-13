
from django.db.models.signals import post_migrate
from django.apps import AppConfig


class EmailerConfig(AppConfig):
    name = 'biostar.emailer'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_group, sender=self)
        post_migrate.connect(init_sub, sender=self)


def init_group(sender, **kwargs):
    """
    Initialize with groups.
    """
    pass


def init_sub(sender, **kwargs):
    """
    Initialize subscriptions.
    """
    pass
