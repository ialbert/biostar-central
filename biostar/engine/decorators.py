from functools import wraps

from django.contrib import messages
from django.urls import reverse
from django.shortcuts import redirect
from . import auth
from django.utils.decorators import available_attrs
from . import models


class object_access:
    """
    Wraps a view and checks user access permissions to an instance.
    Redirects to the url on access error.
    """

    def __init__(self, type, access=models.Access.ADMIN_ACCESS, url='', login_required=False):

        # The object that will be checked for permission.
        self.type = type

        # The required access for the url.
        self.access = access

        # The url to redirect to if the access check fails.
        self.url = url

        # Does the access require a logged in user.
        self.login_required = login_required

    def __call__(self, function, *args, **kwargs):
        """
        Decorator used to test if a user has rights to access an instance
        """

        # Pass function attributes to the wrapper
        @wraps(function, assigned=available_attrs(function))
        def _wrapped_view(request, *args, **kwargs):

            # Each wrapped view must take a numerical id as parameter.
            id = kwargs['id']

            # The user is set in the request.
            user = request.user

            # Fetches the object that will be checked for permissions.
            instance = self.type.objects.filter(pk=id).first()

            # Check for access to the object.
            allow_access = auth.check_obj_access(user=user, instance=instance, request=request, access=self.access, login_required=self.login_required)

            # Access check did not pass, redirect.
            if not allow_access:

                # If there is a redirect url build with the id.
                if self.url:
                    target = reverse(self.url, kwargs=dict(id=id))
                else:
                    target = reverse('project_list')

                return redirect(target)

            # Return the wrapped function.
            return function(request, *args, **kwargs)

        return _wrapped_view


