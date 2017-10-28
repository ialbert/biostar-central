from django.db.models.signals import post_migrate
from django.apps import AppConfig
from . import signals

class EngineConfig(AppConfig):
    name = 'biostar.engine'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(signals.init_site, sender=self)
        post_migrate.connect(signals.init_users, sender=self)
        post_migrate.connect(signals.init_proj, sender=self)



