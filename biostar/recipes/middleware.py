from django.contrib.auth import logout
from django.contrib import messages
from django.conf import settings
from biostar.accounts.models import Profile
from biostar.recipes.auth import detect_cores


#def redirect_middleware():
#    return


def recipes_middleware(get_response):

    def middleware(request):

        user = request.user

        # Banned and suspended users are not allowed
        if user.is_authenticated and user.profile.state in (Profile.BANNED, Profile.SUSPENDED, Profile.SPAMMER):
            messages.error(request, f"Account is {user.profile.get_state_display()}")
            logout(request)

        response = get_response(request)

        origin = detect_cores(request)
        if origin:
            response["Access-Control-Allow-Origin"] = origin

        return response

    return middleware


