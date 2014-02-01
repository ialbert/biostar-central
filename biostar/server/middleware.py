__author__ = 'ialbert'
from django.contrib import messages
from django.contrib.auth import logout

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

        # Check suspended status for users.
        if user.is_authenticated() and user.is_suspended:
            logout(request)
            messages.error(request, 'Sorry, this account has been suspended. Please contact the administrators.')
            return