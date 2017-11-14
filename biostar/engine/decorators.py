from .models import Project
from django.contrib import messages
from django.urls import reverse
from django.shortcuts import redirect
from . import auth


class object_access:

    def __init__(self, instance, owner_only=False, redirect_url=None):
        self.instance = instance
        self.owner_only = owner_only
        self.redirect_url = redirect_url

    def __call__(self, function, *args, **kwargs):
        """
           Decorator used to test if a user has rights to access an instance
       """

        def wrap(request, *args, **kwargs):
            id, user = kwargs['id'], request.user

            query = self.instance.objects.filter(pk=id).first()

            project, instance, allow_access = auth.check_obj_access(user, query,
                                                           owner_only=self.owner_only)

            # Everything has a project associated with it
            if not project:
                messages.error(request, f"Project with id={id} does not exist.")
                return redirect(reverse("project_list"))

            if allow_access:
                return function(request, *args, **kwargs)
            else:
                #TODO: the redirection is generic and annoying
                messages.error(request, f"Access to {self.instance.__name__}={id} denied.")
                return redirect(reverse("project_list"))

        wrap.__doc__ = function.__doc__
        wrap.__name__ = function.__name__
        return wrap



