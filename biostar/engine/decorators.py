from .models import Project
from django.contrib import messages
from django.urls import reverse
from django.shortcuts import redirect
from . import auth


class object_access:

    def __init__(self, instance, owner_only=False):

        self.instance = instance
        self.owner_only = owner_only

    def __call__(self, function, *args, **kwargs):
        """
           Decorator used to test if a user has rights to access a project
       """

        def wrap(request, *args, **kwargs):
            id, user = kwargs['id'], request.user

            query = self.instance.objects.filter(pk=id).first()

            instance, allow_access = auth.check_obj_access(user, query,
                                                           owner_only=self.owner_only)
            # allow_access
            if allow_access:
                return function(request, *args, **kwargs)
            else:
                messages.error(request, "User not allowed to modify project")
                return redirect(reverse("project_list"))

        wrap.__doc__ = function.__doc__
        wrap.__name__ = function.__name__
        return wrap


class project_access:

    def __init__(self, owner_only=False, group_only=False):

        self.owner_only = owner_only
        self.group_only = group_only

    def __call__(self, function, *args, **kwargs):
        """
           Decorator used to test if a user has rights to access a project
       """

        def wrap(request, *args, **kwargs):
            id, user = kwargs['id'], request.user
            project, allow_access = auth.check_project_access(user, id=id,
                                                              owner_only=self.owner_only,
                                                              group_only=self.group_only)

            if not project:
                messages.error(request, f"Project with id={id} does not exist.")
                return redirect(reverse("project_list"))

            if allow_access:
                return function(request, *args, **kwargs)
            else:
                messages.error(request, "User not allowed to access/modify project")
                return redirect(reverse("project_list"))

        wrap.__doc__ = function.__doc__
        wrap.__name__ = function.__name__
        return wrap
