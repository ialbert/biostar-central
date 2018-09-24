
from django.http import JsonResponse
from functools import wraps, partial
from django.utils.decorators import available_attrs
from django.shortcuts import redirect, reverse
from django.contrib import messages


def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return JsonResponse(payload)


ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')


class ajax_error_wrapper:
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, method):
        self.method = method

    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _ajax_view(request, *args, **kwargs):

            if request.method != self.method:
                return ajax_error(f'{self.method} method must be used.')

            if not request.user.is_authenticated:
                return ajax_error('You must be logged in to do that.')

            return func(request, *args, **kwargs)

        return _ajax_view


class object_exists:

    def __init__(self, klass, url=None):

        self.klass = klass
        self.url = url or "/"

    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _wrapper_function(request, *args, **kwargs):

            # Each wrapped view must take an alphanumeric uid as parameter.
            uid = kwargs.get('uid')

            if not self.klass.objects.get_all(uid=uid).exists():
                messages.error(request, "Object does not exist.")
                return redirect(self.url)

            return func(request, *args, **kwargs)

        return _wrapper_function