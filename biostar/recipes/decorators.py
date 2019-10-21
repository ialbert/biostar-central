from functools import wraps

from django.contrib import messages
from django.shortcuts import redirect, reverse
from django.utils.decorators import available_attrs
from django.conf import settings
from django.http import HttpResponse

from rest_framework import status
from . import models, auth

# Share the logger with models
logger = models.logger


class read_access:
    """
    Controls READ level access to urls.
    """

    def __init__(self, obj_type):
        self.type = obj_type
        self.fallback_url = lambda : reverse("project_list")

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
                return redirect(self.fallback_url())

            # Get project for the instance.
            project = instance.project

            # Allow read access to public projects
            if project.is_public:
                return function(request, *args, **kwargs)

            # Anonymous users may not access non public projects.
            if user.is_anonymous:
                messages.error(request, f"You must be logged in to access object id {uid}")
                return redirect(self.fallback_url())

            # Check the presence of READ or WRITE access
            read_or_write = [models.Access.READ_ACCESS, models.Access.WRITE_ACCESS]
            access = models.Access.objects.filter(user=user, project=project, access__in=read_or_write).first()

            # Project owners may read their project.
            if access or project.owner == user:
                return function(request, *args, **kwargs)

            # Deny access by default.
            msg = auth.access_denied_message(user=user, needed_access=models.Access.READ_ACCESS)
            messages.error(request, msg)
            return redirect(self.fallback_url())

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
        Decorator used to tested if a user has rights to access an instance
        """

        # Pass function attributes to the wrapper
        @wraps(function, assigned=available_attrs(function))
        def _wrapped_view(request, *args, **kwargs):
            # Each wrapped view must take an alphanumeric uid as parameter.
            uid = kwargs.get('uid')
            user = request.user

            # Fetches the object that will be checked for permissions.
            instance = self.type.objects.filter(uid=uid).first()
            if not instance:
                messages.error(request, f"Object id {uid} does not exist.")
                return redirect(reverse("project_list"))

            # Build redirect url
            if self.fallback_view:
                target = reverse(self.fallback_view, kwargs=dict(uid=uid))
            else:
                target = request.GET.get("next") or instance.url()

            # The project that corresponds to the instance.
            project = instance.project

            # Check write to an object.
            access = auth.has_write_access(user=user, project=project)

            # Project owners may write their project.
            if access:
                return function(request, *args, **kwargs)

            msg = auth.access_denied_message(user=user, needed_access=models.Access.WRITE_ACCESS)
            messages.error(request, msg)
            return redirect(target)

        return _wrapped_view


def parse_api_key(request):

    empty = ""
    if request.method == "PUT":
        return request.data.get("k", empty)
    elif request.method == "GET":
        return request.GET.get("k", empty)

    return empty


class require_api_key:
    """
    Require the api key when private project or PUT request for an api view
    """

    def __init__(self, type):
        self.type = type

    def __call__(self, func, *args, **kwargs):

        # Pass function attributes to the wrapper
        @wraps(func, assigned=available_attrs(func))
        def _api_view(request, *args, **kwargs):
            # Each wrapped view must take an alphanumeric uid as parameter.
            api_key = parse_api_key(request=request)

            # Get the Recipe uid
            uid = kwargs.get("uid")
            obj = self.type.objects.filter(uid=uid).first()

            if not obj:
                msg = dict(error="Object does not exist.")
                return HttpResponse(content=msg, status=status.HTTP_404_NOT_FOUND)

            # All PUT requests will require an API key.
            if request.method == "PUT" and settings.API_KEY != api_key:
                msg = dict(error="API key is required for all PUT requests.")
                return HttpResponse(content=msg, status=status.HTTP_401_UNAUTHORIZED)

            # API key required when asking for GET requests of private recipes.
            elif request.method == "GET" and settings.API_KEY != api_key and obj.project.is_private:
                msg = dict(error="Private recipes can not be accessed without an API key param (?k=).")
                return HttpResponse(content=msg, status=status.HTTP_401_UNAUTHORIZED)

            return func(request, *args, **kwargs)

        return _api_view

