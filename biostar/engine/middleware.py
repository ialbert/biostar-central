from django.contrib import auth
from django.contrib import messages

from biostar.accounts.models import Profile


def engine_middleware(get_response):

    def middleware(request):

        user = request.user

        # Banned and suspended users are not allowed
        if user.is_authenticated and user.profile.state in (Profile.BANNED, Profile.SUSPENDED):
            messages.error(request, f"Account is {user.profile.get_state_display()}")
            auth.logout(request)

        response = get_response(request)
        # Can process response here after its been handled by the view

        return response

    return middleware

