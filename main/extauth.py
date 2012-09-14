# external authorization
import random
from django.conf import settings
from main.server import models
from django.contrib import messages

# this is the external authenticator
from django.contrib.auth import authenticate, login, logout

def dummy_auth(request):
    num   = random.randint(1,3)
    name  = "John Doe %s" % num
    email = "john%s@email.com" % num
    value = dict(name=name, email=email)
    return value        

# this will be the function that will need to be called to authenticate
external_authenticator_func = globals()[settings.EXTERNAL_AUTHENTICATOR_FUNC]

class ExternalAuthenticator(object):
    """
    Attempts to authenticate a user from an external source
    """
    def process_request(self, request):
        user = request.user
        anon = not user.is_authenticated()
        info = external_authenticator_func(request)
        if anon and info:
            # attempt to authenticate here
            name, email = info['name'], info['email']
            user = models.User.objects.filter(email=email)
            if not user:
                username = models.make_uuid()
                user = models.User(email=email, username=username)
                user.save()
                user.profile.display_name = name
                user.profile.save()
            else:
                assert len(user) == 1, "duplicated emails, why?"
                user = user[0]
            password = models.make_uuid()
            user.set_password(password)
            user.save()
            user = authenticate(username=user.username, password=password)
            messages.info(request, "Login complete.")
            login(request=request, user=user)
            