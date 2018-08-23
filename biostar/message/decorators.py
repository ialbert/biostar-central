from functools import wraps
from django.utils.decorators import available_attrs
from django.shortcuts import redirect, reverse
from django.db.models import Q
from django.contrib import messages

from .const import *
from . import models


class message_access:

    def __init__(self, access_to=INBOX):

        self.access_to = access_to

    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _wrapped_view(request, *args, **kwargs):

            # Each wrapped view must take an alphanumeric uid as parameter.
            uid = kwargs.get('uid')
            user = request.user
            if user.is_anonymous:
                messages.error(request, "You need to be logged in to view messages.")
                return redirect("/")

            # Only recipients get access to the inbox and sender to outbox.
            if self.access_to == OUTBOX:
                msg_user = Q(sender=user)
            else:
                msg_user = Q(recipient=user)

            # User needs to be the sender or recipient to view message
            message = models.Message.objects.filter(msg_user, uid=uid).exists()

            if message:
                return func(request, *args, **kwargs)

            messages.error(request, "Object does not exist")
            return redirect(reverse("inbox"))

        return _wrapped_view

