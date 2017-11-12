from .models import Project
from django.contrib import messages
from django.urls import reverse
from django.shortcuts import redirect


def access_project(function):
    """
    Decorator used to test if a user has rights to access a project
    """
    def wrap(request, *args, **kwargs):

        project = Project.objects.filter(pk=kwargs['id']).first()

        if project.owner == request.user or request.user.is_superuser:
            return function(request, *args, **kwargs)
        else:
            messages.error(request, "User not allowed to modify project")
            return redirect(reverse("project_view", kwargs=dict(id=project.id)))

    wrap.__doc__ = function.__doc__
    wrap.__name__ = function.__name__
    return wrap

