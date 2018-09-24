from biostar.utils.decorators import *
from functools import wraps
from django.utils.decorators import available_attrs
from django.shortcuts import redirect, reverse
from django.contrib import messages
from . import const


def protect_private_topics(func):

    @wraps(func, assigned=available_attrs(func))
    def _wrapped_view(request, *args, **kwargs):

        topic = kwargs.get("topic", request.GET.get("active", const.LATEST))

        is_private_topic = topic in const.PRIVATE_TOPICS

        if request.user.is_anonymous and is_private_topic:
            messages.error(request, f"You must be logged in to view that topic.")
            return redirect(reverse("post_list"))

        return func(request, *args, **kwargs)

    return _wrapped_view


