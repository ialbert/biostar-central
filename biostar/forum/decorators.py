from functools import wraps, partial
from django.utils.decorators import available_attrs
from django.shortcuts import redirect, reverse
from django.contrib import messages
from django.http import JsonResponse

from . import const


class object_exists:

    def __init__(self, klass, url=None):

        self.klass = klass
        self.url = url or "post_list"

    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _wrapper_function(request, *args, **kwargs):

            # Each wrapped view must take an alphanumeric uid as parameter.
            uid = kwargs.get('uid')

            if not self.klass.objects.get_all(uid=uid).exists():
                messages.error(request, "Object does not exist.")
                return redirect(reverse(self.url))

            return func(request, *args, **kwargs)

        return _wrapper_function


def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return JsonResponse(payload)


ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')


class ajax_error_wrapper(object):
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, func):
        self.func = func

    def __call__(self, request):
        try:
            if request.method != 'POST':
                return ajax_error('POST method must be used.')

            if not request.user.is_authenticated:
                return ajax_error('You must be logged in to do that')

            return self.func(request)

        except Exception as exc:
            return ajax_error('Error: %s' % exc)


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


