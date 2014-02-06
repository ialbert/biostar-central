__author__ = 'ialbert'

def user_permissions(request, target):
    """
    Sets permission attributes on a user
    """
    user = request.user
    is_editable = has_ownership = False

    if user.is_authenticated():

        if user == target:
            has_ownership = is_editable = True

        elif user.is_administrator and not target.is_administrator:
            is_editable = True

        elif user.is_moderator and not target.is_moderator:
            is_editable = True

    target.is_editable = is_editable
    target.has_ownership = has_ownership

    return target