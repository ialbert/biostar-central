from functools import wraps
from django.utils.decorators import available_attrs
from django.shortcuts import redirect, reverse
from django.contrib import messages




class object_exists:

    def __init__(self, klass, url=None):

        self.klass = klass
        self.url = url or "post_list"


    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _wrapper_function(request, *args, **kwargs):

            # Each wrapped view must take an alphanumeric uid as parameter.
            uid = kwargs.get('uid')

            if not self.klass.objects.filter(uid=uid).exists():
                messages.error(request, "Object does not exist.")
                return redirect(reverse(self.url))

            return func(request, *args, **kwargs)

        return _wrapper_function
