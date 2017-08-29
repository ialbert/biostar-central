from django.contrib.auth import get_user_model
from django.contrib.auth.backends import ModelBackend


class DualLoginModelBackend(ModelBackend):
    
    def authenticate(self, username=None, password=None):
        usermodel = get_user_model()

        if '@' in username:
            kwargs = {'email': username}
        else:
            kwargs = {'username': username}

        try:
            user = usermodel.objects.get(**kwargs)
            if user.check_password(password):
                return user
        except usermodel.DoesNotExist:
            return None

    def get_user(self, user_id):
        usermodel = get_user_model()
        try:
            return usermodel.objects.get(pk=user_id)
        except usermodel.DoesNotExist:
            return None




