from functools import wraps

from django.urls import reverse
from django.shortcuts import redirect
from django.utils.decorators import available_attrs
from biostar.accounts.models import Profile
from . import auth

from . import models



class object_access:
    """
    Wraps a view and checks user access permissions to an instance.
    Redirects to the url on access error.
    """


    def __init__(self, type, access=models.Access.WRITE_ACCESS, url='', login_required=False):

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

            # Each wrapped view must take an alphanumeric uid as parameter.
            uid = kwargs.get('uid')

            # The user is set in the request.
            user = request.user

            # Fetches the object that will be checked for permissions.
            instance = self.type.objects.filter(uid=uid).first()

            # Check for access to the object.
            allow_access = auth.check_obj_access(user=user, instance=instance, request=request, access=self.access,
                                                 login_required=self.login_required)
            # Access check did not pass, redirect.
            if not allow_access:

                # If there is a redirect url build with the uid.
                if self.url and uid:
                    target = reverse(self.url, kwargs=dict(uid=uid))
                else:
                    target = reverse('project_list')

                return redirect(target)

            # Return the wrapped function.
            return function(request, *args, **kwargs)

        return _wrapped_view


