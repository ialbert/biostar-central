from django.contrib.syndication.views import Feed
from django.shortcuts import get_object_or_404
from main.server import models, const, html

class LatestEntriesFeed(Feed):
    title = "Biostars.org latest"
    link = "/"
    description = "Latest 25 posts from the Biostar server"

    def items(self):
        return models.Post.objects.filter(type=const.POST_QUESTION).order_by('-creation_date')[:25]

    def item_title(self, item):
        return item.title

    def item_description(self, item):
        return item.html[:1000]

class NotificationFeed(Feed):
    title = "Biostar notifications"
    link = "/"
    description = "Latest 25 notification for a given user"

    def title(self, obj):
        return " Notifications for %s" % obj.profile.display_name

    def get_object(self, request, uuid):
        obj = get_object_or_404(models.User, profile__uuid=uuid)
        return obj
        
    def items(self, obj):
        return models.Note.objects.filter(target=obj).select_related('sender', 'sender__profile').order_by('-date')[:25]

    def item_title(self, item):
        return "From %s" % item.sender.profile.display_name

    def item_description(self, item):
        return item.html[:1000]

class MyTagsFeed(Feed):
    title = "Biostar MyTags"
    link = "/"
    description = "Latest posts matching your tags"

    def title(self, obj):
        return "Post matching tags for %s" % obj.profile.display_name

    def get_object(self, request, uuid):
        obj = get_object_or_404(models.User, profile__uuid=uuid)
        return obj
        
    def items(self, obj):
        tags  = obj.profile.my_tags.split(' ')
        posts = models.query_by_tags(user=obj, tags=tags).order_by('-creation_date')
        return posts[:15]

    def item_title(self, item):
        return item.title

    def item_description(self, item):
        return item.content[:1000]