import random

from django.apps import AppConfig
from django.conf import settings

from django.db.models.signals import post_migrate


class EmailerConfig(AppConfig):
    name = 'biostar.emailer'

    def ready(self):
        # Triggered upon app initialization.
        pass
