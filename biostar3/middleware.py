from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib import messages
from django.conf import settings
from django.contrib.auth import get_user_model

User = get_user_model()

def add_user_attributes(user):
    """
    Ensures that anyonymous users have attributes that autheniticated users have.
    This greatly simplfies templating. Will mutate data for anonymous users.
    """
    if not user.is_authenticated():
        user.is_moderator = user.is_admin = False

class GlobalMiddleware(object):
    """Performs tasks that are applied on every request"""

    def process_request(self, request):

        # Ensure that anonymous users have all required attributes.
        add_user_attributes(request.user)
