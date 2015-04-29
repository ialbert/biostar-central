from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib import admin
from .models import User, Post, Blog, BlogPost, FlatPage

admin.site.register(User)

admin.site.register(Post)

admin.site.register(Blog)

admin.site.register(BlogPost)

admin.site.register(FlatPage)

