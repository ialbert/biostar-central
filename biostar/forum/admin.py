
from django.contrib import admin
from biostar.forum.models import Post, Subscription, Vote, Log, SharedLink


@admin.register(Post)
class PostAdmin(admin.ModelAdmin):
    list_display = ('uid', 'title', 'type', 'author', 'lastedit_date')
    ordering = ['type', 'rank']
    fieldsets = (
        (None, {'fields': ('title',)}),
        ('Attributes', {'fields': ('type', 'status', 'author', 'sticky',)}),
        ('Content', {'fields': ('content', )}),
    )
    search_fields = ('title', 'author__profile__name', 'uid')


@admin.register(Vote)
class VoteAdmin(admin.ModelAdmin):
    list_display = ('author', 'post', 'type', 'date')
    ordering = ['-date']
    search_fields = ('post__title', 'author__profile__name', 'uid')


@admin.register(Subscription)
class SubscriptionAdmin(admin.ModelAdmin):
    search_fields = ('user__profile__name', 'user__email', 'uid')
    list_select_related = ["user", "post"]

@admin.register(Log)
class LoggerAdmin(admin.ModelAdmin):
    list_display = ('action', 'user', 'date', 'text', "ipaddr")
    ordering = ['-date']

    fieldsets = (("Details",
                  {"fields": ("user", 'action', 'text')}
                  ),)

    search_fields = ('user__email', 'user__profile__name', 'user__profile__uid', 'log_text')


@admin.register(SharedLink)
class SharedLinkAdmin(admin.ModelAdmin):
    list_display = ('url',  'status', 'author', 'title')
    search_fields = ('author__email', 'author__profile__name', 'url', 'title')

    pass
