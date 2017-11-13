from django.contrib.auth.models import User
from django.contrib import messages
from django.urls import reverse
from django.shortcuts import redirect

def user_level_access(function):
    """
       Decorator used to test if user being queried is the same as request.user
   """

    def wrap(request, *args, **kwargs):
        id = kwargs['id']
        user = User.objects.filter(pk=id).first()

        if (user == request.user):
            return function(request, *args, **kwargs)
        elif request.user.is_authenticated:
            messages.error(request, "User not allowed to access profile")
            return redirect(reverse("profile", kwargs=dict(id=user.id)))
        else:
            messages.error(request, "User not allowed to access profile")

    wrap.__doc__ = function.__doc__
    wrap.__name__ = function.__name__
    return wrap