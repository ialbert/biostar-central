import os
from django.db.models.signals import post_migrate
from django.apps import AppConfig
from django.conf import settings


class PlanetConfig(AppConfig):
    name = 'biostar.planet'

    def ready(self):
        #from . import signals
        # Triggered upon app initialization.
        #post_migrate.connect(init_awards, sender=self)
        #post_migrate.connect(init_digest, sender=self)
        #post_migrate.connect(init_planet, sender=self)
        pass