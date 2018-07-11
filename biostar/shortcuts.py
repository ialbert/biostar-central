from django.urls import reverse as main_reverse
from django.contrib import messages


def reverse(view, request=None, **kwargs):

    try:
        url = main_reverse(view, kwargs=kwargs)
    except Exception as exc:
        url = "/"
        if request:
            messages.error(request, f"Error reversing: {view}")

    return url





