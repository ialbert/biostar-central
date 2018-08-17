from functools import wraps
from django.utils.decorators import available_attrs
from django.shortcuts import redirect, reverse
from django.db.models import Q


from . import models


def message_access(func):

    @wraps(func, assigned=available_attrs(func))
    def _wrapped_view(request, *args, **kwargs):

        # Each wrapped view must take an alphanumeric uid as parameter.
        uid = kwargs.get('uid')
        user = request.user

        # User needs to be the sender or recipient to view message
        message = models.Message.objects.filter(Q(sender=user) | Q(recipient=user), uid=uid).exists()

        # Only sender and recipient see the message.
        if message:
            return func(request, *args, **kwargs)

        return redirect(reverse("message_list"))

    return _wrapped_view

