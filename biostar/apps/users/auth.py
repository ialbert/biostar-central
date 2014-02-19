__author__ = 'ialbert'

def user_permissions(request, target):
    """
    Sets permission attributes on a user.

    The is_staff attribute is internal to Django. It allows a user
    to log into the Django admin. The role should be reserved to only
    those that manage the servers.
    """
    user = request.user

    has_ownership = is_editable = False

    if not user.is_authenticated():
        # Anonymous users cannot do anything.
        has_ownership = target.is_editable = False
    elif user.is_staff:
        # Django level staff have full access.
        has_ownership = is_editable = True
    elif user == target:
        # A user has full access to themselves.
        has_ownership = is_editable = True
    elif target.is_staff or target.is_administrator:
        # Cannot edit admins or staff.
        is_editable = False
    elif user.is_moderator:
        # User has to be a moderator.
        is_editable = True

    # Apply the attributes
    target.has_ownership = has_ownership
    target.is_editable = is_editable

    return target