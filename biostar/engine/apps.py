from django.apps import AppConfig

class EngineConfig(AppConfig):
    name = 'biostar.engine'

    def ready(self):
        # Triggered upon app initialization.
        pass



