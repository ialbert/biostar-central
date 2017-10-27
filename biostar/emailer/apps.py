
from django.db.models.signals import post_migrate
from django.apps import AppConfig
from . import signals



class EmailerConfig(AppConfig):
    name = 'emailer'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(signals.init_site, sender=self)
