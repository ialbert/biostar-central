from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib import messages
from django.conf import settings
from django.contrib.auth import get_user_model

User = get_user_model()

def modify_request(request):
    """
    Each request will be altered to contain settings that the site expects.

    1. Ensures that anyonymous users have attributes that autheniticated users have.
    This greatly simplfies templating. Will mutate data for anonymous users.

    2. Sets the subdomain of the current request. It is used to filter content by.

    """
    user = request.user

    if not user.is_authenticated():
        user.is_moderator = user.is_admin = False

    domain = request.META['HTTP_HOST']
    request.subdomain = domain.split('.')[0]

class GlobalMiddleware(object):
    """Performs tasks that are applied on every request"""

    def process_request(self, request):

        # Ensures that requests have all the information needed.
        modify_request(request)

