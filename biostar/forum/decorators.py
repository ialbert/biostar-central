
from django.shortcuts import redirect, reverse
from django.contrib import messages

from .models import Post



def object_exists(func):

    def _wrapper_function(request, uid):

        if not Post.objects.filter(uid=uid).exists():
            messages.error(request, "Object does not exist.")
            return redirect(reverse("post_list"))

        return func(request, uid)


    return _wrapper_function