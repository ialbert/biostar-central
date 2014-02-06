__author__ = 'ialbert'

def user_permissions(request, target):
    """
    Sets permission attributes on a user.

    The is_staff attribute is internal to Django. It allows a user
    to log into the Django admin. The role should be reserved to only
    those that manage the servers.
    """
    user = request.user
    is_editable = has_ownership = False

    if user.is_authenticated():

        if user == target or user.is_staff:
            has_ownership = is_editable = True

        elif user.is_administrator and (not target.is_administrator and not target.is_staff):
            is_editable = True

        elif user.is_moderator and (not target.is_moderator and not target.is_staff):
            is_editable = True

    target.is_editable = is_editable
    target.has_ownership = has_ownership

    return target