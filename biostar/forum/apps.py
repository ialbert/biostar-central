import logging
from django.apps import AppConfig

logger = logging.getLogger('engine')


class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        from . import signals
        # Triggered upon app initialization.