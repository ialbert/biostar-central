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
            instance = self.type.objects.get_all(uid=uid).first()

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
            messages.error(request, f"Read access denied to object id: {uid}")
            return redirect(self.fallback_url)

        return wrapper


class write_access:
    """
    Controls WRITE level access to urls.
    """

    def __init__(self, type, fallback_view=None):
        self.type = type
        self.fallback_view = fallback_view

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
            if user.is_anonymous:
                messages.error(request, "You need to be logged in to preform actions.")
                return redirect(reverse("project_list"))

            # Fetches the object that will be checked for permissions.
            instance = self.type.objects.get_all(uid=uid).first()
            if not instance:
                messages.error(request, f"Object id {uid} does not exist.")
                return redirect(reverse("project_list"))

            project = instance.project
            access = models.Access.objects.filter(user=user, project=project, access=models.Access.WRITE_ACCESS).first()

            # Project owners may write their project.
            if access or instance.project.owner == user:
                return function(request, *args, **kwargs)

            # Build redirect url
            if self.fallback_view:
                target = reverse(self.fallback_view, kwargs=dict(uid=uid))
            else:
                target = request.GET.get("next") or instance.url()

            messages.error(request, f"Write access denied to object id: {uid}")
            return redirect(target)

        return _wrapped_view
