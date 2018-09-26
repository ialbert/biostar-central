from functools import wraps
from django.contrib import messages

from django.shortcuts import redirect
from django.utils.decorators import available_attrs

from biostar.utils.shortcuts import reverse
from . import auth

from . import models

# Share the logger with models
logger = models.logger

class read_access:
    """
    Controls READ level access to urls.
    """

    def __init__(self, type):
        self.type = type
        self.fallback_url = reverse("project_list")

    def __call__(self, function, *args, **kwargs):
        # Pass function attributes to the wrapper
        @wraps(function, assigned=available_attrs(function))
        def wrapper(request, *args, **kwargs):

            # Each wrapped view must take an alphanumeric uid as parameter.
            uid = kwargs.get('uid')

            # The user is set in the request.
            user = request.user

            # Fetches the object that will be checked for permissions.
            instance = self.type.objects.filter(uid=uid).first()

            # Object does not exist.
            if not instance:
                messages.error(request, f"Object id {uid} does not exist")
                return redirect(self.fallback_url)

            # Get project for the instance.
            project = instance.project

            # Allow read access to public projects
            if project.is_public:
                return function(request, *args, **kwargs)

            # Anonymous users may not access non public projects.
            if user.is_anonymous:
                messages.error(request, f"You must be logged in to access object id {uid}")
                return redirect(self.fallback_url)

            # Check the presence of READ or WRITE access
            read_or_write = [models.Access.READ_ACCESS, models.Access.WRITE_ACCESS]
            access = models.Access.objects.filter(user=user, project=project, access__in=read_or_write).first()

            # Project owners may read their project.
            if access or project.owner == user:
                return function(request, *args, **kwargs)

            # Deny access by default.
            messages.error(request, f"Read access denied to object id {uid}")
            return redirect(self.fallback_url)

        return wrapper


class write_access:
    """
    Controls WRITE level access to urls.
    """

    def __init__(self, type):
        self.type = type
        self.fallback_url = reverse("project_list")

    def __call__(self, function, *args, **kwargs):
        """
        Decorator used to test if a user has rights to access an instance
        """

        # Pass function attributes to the wrapper
        @wraps(function, assigned=available_attrs(function))
        def _wrapped_view(request, *args, **kwargs):

            pass


class owner_only:
    """
    Controls privileges left to owner
    """
    def __init__(self, type):
        self.type = type
        self.fallback_url = reverse("project_list")

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

            # Project owners may read their project.
            if instance.owner == user or instance.project.owner == user:
                return function(request, *args, **kwargs)

            messages.error(request, "You have to be the owner to preform this action.")
            return redirect(self.fallback_url)

        return _wrapped_view


class object_access:
    """
    Wraps a view and checks user access permissions to an instance.
    Redirects to the url on access error.
    """

    def __init__(self, type, access=models.Access.WRITE_ACCESS, url='', role=None,
                 login_required=False, show_deleted=True):

        # The object that will be checked for permission.
        self.type = type

        # The required access for the url.
        self.access = access

        # The url to redirect to if the access check fails.
        self.url = url

        # Does the access require a logged in user.
        self.login_required = login_required

        # Required role of user. Given precedent over access
        self.role = role

        # Allow to show deleted objects in the view
        self.show_deleted = show_deleted

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
            instance = self.type.objects.get_all(uid=uid).first()

            # Check for access to the object.
            allow_access = auth.check_obj_access(user=user, instance=instance, request=request, access=self.access,
                                                 login_required=self.login_required, role=self.role)
            if self.url and uid:
                target = reverse(self.url, request=request, kwargs=dict(uid=uid))
            else:
                target = reverse('project_list', request=request)

            target = request.GET.get("next") or target
            # Access check did not pass, redirect.
            if not allow_access:
                return redirect(target)
            elif not self.show_deleted and instance.deleted:
                messages.error(request, "Can not preform this action while object is in the recycle bin.")
                return redirect(target)

            # Return the wrapped function.
            return function(request, *args, **kwargs)

        return _wrapped_view
