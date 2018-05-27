from django.apps import AppConfig
from django.conf import settings



class ForumConfig(AppConfig):
    name = 'biostar.forum'

    def ready(self):
        # Triggered upon app initialization.
        pass




def init_post(sender,  **kwargs):



    from django.contrib.auth import get_user_model


    User = get_user_model()

    name, email = settings.ADMINS[0]

    user = User.objects.filter(email=email).first()








    return