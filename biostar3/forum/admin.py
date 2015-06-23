from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib import admin
from . import models

class ProfileInline(admin.StackedInline):
    model = models.Profile
    fields = [ "info", "scholar", "twitter_id", "website" ]

class FlatPageInline(admin.StackedInline):
    model = models.FlatPage


@admin.register(models.User)
class UserAdmin(admin.ModelAdmin):
    list_select_related = [ 'profile' ]
    search_fields = ['email', 'name']
    inlines = [
     ProfileInline
    ]
    fields = ('email', 'name', 'status')

@admin.register(models.Post)
class PostAdmin(admin.ModelAdmin):
    list_select_related = [ 'author' ]
    search_fields = ['title', 'id']

    inlines = [
        FlatPageInline
    ]

    #fields = ('email', 'name', 'status')


admin.site.register(models.Blog)

admin.site.register(models.BlogPost)

admin.site.register(models.FlatPage)
