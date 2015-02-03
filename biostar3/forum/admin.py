from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib import admin
from .models import User, Post

admin.site.register(User)

admin.site.register(Post)
