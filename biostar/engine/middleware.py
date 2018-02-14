from django.utils.deprecation import MiddlewareMixin
from biostar.accounts.models import Profile
from django.http import HttpResponseRedirect

from django.contrib import messages

import os


class EngineMiddleware(MiddlewareMixin):

    def process_request(self, request):
        "Redirect to root url if the user is not trusted"

        user = request.user
        if user.is_anonymous:
            # Continue toProcess other middleware
            return None

        is_root_url = request.get_full_path() == "/"

        # Request only allowed if the user is 'new' or 'trusted'
        if user.profile.state in (Profile.NEW, Profile.TRUSTED) or is_root_url :
            return None

        messages.error(request, f"Account is {user.profile.get_state_display()}")
        return HttpResponseRedirect("/")


