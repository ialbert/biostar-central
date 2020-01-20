from django.contrib import admin
from biostar.planet.models import Blog


@admin.register(Blog)
class BlogAdmin(admin.ModelAdmin):
    list_display = ('title', 'feed', 'list_order', 'link')

    ordering = ['list_order', 'active']

    fieldsets = (
        (None, {'fields': ('title', 'desc', 'feed', 'link', 'active', 'list_order')}),
    )

    search_fields = ('title', 'desc', 'link', 'feed', 'list_order')
