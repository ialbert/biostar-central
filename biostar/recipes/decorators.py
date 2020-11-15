from functools import wraps

from django.contrib import messages
from django.shortcuts import redirect, reverse
from django.conf import settings
from django.http import HttpResponse
from django.http import QueryDict

from biostar.accounts.models import User
from . import models, auth

# Share the logger with models
logger = models.logger


class read_access:
    """
    Controls READ level access to urls.
    """

    def __init__(self, type, allowed_cors=None, strict=False, fallback_view="", login_required=False, json=False):
        self.type = type
        self.allowed_cors = allowed_cors
        self.login_required = login_required
        self.fallback_view = fallback_view
        # Strict policy enforced so public projects still get their access checked.
        self.strict = strict

    def __call__(self, function, *args, **kwargs):
        # Pass function attributes to the wrapper
        @wraps(function)
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
                return redirect("project_list")

            # Get project for the instance.
            project = instance.project

            # Build redirect url
            if self.fallback_view:
                target = reverse(self.fallback_view, kwargs=dict(uid=uid))
            else:
                target = request.GET.get("next") or "project_list"

            # Allow read access for allowed CORS websites.
            origin = auth.detect_cores(request)
            if self.allowed_cors and self.allowed_cors in origin:
                return function(request, *args, **kwargs)

            # Authenticated users required for to view this page.
            if self.login_required and request.user.is_anonymous:
                messages.error(request, "You need to be logged in.")
                return redirect(target)

            # Allow read access to public projects
            if project.is_public:
                return function(request, *args, **kwargs)

            # Anonymous users may not access non public projects.
            if user.is_anonymous:
                messages.error(request, f"You must be logged in to access object id {uid}")
                return redirect("project_list")

            # Check the presence of READ or WRITE access
            readable = auth.is_readable(user=user, obj=project, strict=self.strict)

            # Project owners may read their project.
            if readable or project.owner == user:
                return function(request, *args, **kwargs)

            # Deny access by default.
            msg = auth.access_denied_message(user=user, needed_access=models.Access.READ_ACCESS)
            messages.error(request, msg)
            return redirect("project_list")

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
        @wraps(function)
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
            access = auth.is_writable(user=user, project=project)

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
        return QueryDict(request.body)
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
        @wraps(func)
        def _api_view(request, *args, **kwargs):
            # Each wrapped view must take an alphanumeric uid as parameter.
            api_key = parse_api_key(request=request)

            # Get the Recipe uid
            uid = kwargs.get("uid")
            obj = self.type.objects.filter(uid=uid).first()

            if not obj:
                msg = dict(error="Object does not exist.")
                return HttpResponse(content=msg, status=404)

            # All PUT requests will require an API key.
            if request.method == "PUT" and settings.API_KEY != api_key:
                msg = dict(error="API key is required for all PUT requests.")
                return HttpResponse(content=msg, status=401)

            # API key required when asking for GET requests of private recipes.
            elif request.method == "GET" and settings.API_KEY != api_key and obj.project.is_private:
                msg = dict(error="Private recipes can not be accessed without an API key param (?k=).")
                return HttpResponse(content=msg, status=401)

            return func(request, *args, **kwargs)

        return _api_view


def token_access(klass, allow_create=False):
    """
    Check users access to an object using their token.
    """

    def __outer__(func):

        @wraps(func)
        def __wrapper__(request, *args, **kwargs):

            # Get the token from the request data
            token = auth.get_token(request=request)
            # Find the target user.
            user = User.objects.filter(profile__token=token).first()

            # Get the object user is trying to access
            uid = request.GET.get('uid', request.POST.get('uid', ''))
            obj = klass.objects.filter(uid=uid).first()

            if not obj:
                # Allow users to create.
                if allow_create and user:
                    return func(request, *args, **kwargs)
                return HttpResponse(content="Object does not exist.", status=404)

            project = obj.project

            # User token required for private projects.
            if not user and project.is_private:
                return HttpResponse(content="Token does not belong to any user.", status=404)

            # GET requests require read access
            if request.method == "GET":
                acc = auth.is_readable(user=user, obj=obj, strict=False)
            # PUT and POST requests require write access
            else:
                acc = auth.is_writable(user=user, project=project)

            # User passes test.
            if acc:
                return func(request, *args, **kwargs)

            return HttpResponse(content="User does not have access to preform that action.", status=404)

        return __wrapper__

    return __outer__
