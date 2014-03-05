__author__ = 'ialbert'
from django.conf import settings
from datetime import datetime, timedelta

def user_permissions(request, target):
    """
    Sets permission attributes on a user.

    The is_staff attribute is internal to Django. It allows a user
    to log into the Django admin. The role should be reserved to only
    those that manage the servers.
    """
    user = request.user

    # The user is the target.
    has_ownership = is_editable = False

    if user.is_authenticated():

        if user == target:
            # A user can do anything on their own account.
            has_ownership = is_editable = (user == target)

        elif target.is_administrator:
            # Admins cannot be moderated.
            is_editable = False

        elif user.is_administrator:
            # Admins can edit other users
            is_editable = True

        elif user.is_moderator and not target.is_moderator:
            # A moderator can edit other non-moderators.
            is_editable = True

    # Apply the attributes
    target.has_ownership = has_ownership
    target.is_editable = is_editable

    return target