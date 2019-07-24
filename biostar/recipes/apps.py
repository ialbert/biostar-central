from django.apps import AppConfig

class EngineConfig(AppConfig):
    name = 'biostar.recipes'

    def ready(self):
        # Triggered upon app initialization.
        pass



