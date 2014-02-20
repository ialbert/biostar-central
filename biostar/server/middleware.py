__author__ = 'ialbert'
from django.contrib import messages
from django.contrib.auth import logout
from django.conf import settings
import hmac, logging
from django.contrib.auth import authenticate, login, logout
from biostar.apps.users.models import User

logger = logging.getLogger(__name__)

def valid_external_login(request):
    "Attempts to perform an external login"

    for name, key in settings.EXTERNAL_AUTH:
        value = request.COOKIES.get(name)
        if value:
            try:
                email, digest1 = value.split(":")
                digest2 = hmac.new(key, email).hexdigest()
                if digest1 != digest2:
                    raise Exception("digests do not match: %s vs %s" % (digest1, digest2))
            except Exception, exc:
                logger.error(exc)
                return None

            # If we made it this far the data is valid.
            password = settings.SECRET_KEY + email

            user, flag = User.objects.get_or_create(email=email)
            if flag:
                logger.info("created user %s" % user.email)

            # Regular login is not allowed.
            user.set_password(password)
            user.save()

            # Authenticate.
            user = authenticate(username=user.email, password=password)
            login(request=request, user=user)
            return True

    return False

class Visit(object):
    """
    Sets visit specific parameters on objects.
    """

    def process_request(self, request):


        user = request.user

        # Add attributes to anonymous users.
        if not user.is_authenticated():

            # This attribute is required inside templates.
            user.is_moderator = False

            if settings.EXTERNAL_AUTH and valid_external_login(request):
                messages.success(request, "Login completed")

        # Check suspended status for users.
        if user.is_authenticated() and user.is_suspended:
            logout(request)
            messages.error(request, 'Sorry, this account has been suspended. Please contact the administrators.')
            return