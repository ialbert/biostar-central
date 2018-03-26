from django.apps import AppConfig

class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        # Triggered upon app initialization.
        pass