from django.db.models.signals import post_migrate
from django.apps import AppConfig
from .signals import init_proj, init_site, init_users, init_groups


class EngineConfig(AppConfig):
    name = 'engine'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_site, sender=self)
        post_migrate.connect(init_groups, sender=self)
        post_migrate.connect(init_users, sender=self)
        post_migrate.connect(init_proj, sender=self)




