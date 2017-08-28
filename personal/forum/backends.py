from django.contrib.auth import get_user_model
from django.contrib.auth.backends import ModelBackend

User = get_user_model()


class EmailModelBackend(ModelBackend):
   

    def authenticate(self, username=None, password=None):
         return NotImplemented 

    def get_user(self, user_id=None):
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None
