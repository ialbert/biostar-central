from django.contrib.syndication.views import Feed
from django.shortcuts import get_object_or_404
from main.server import models, const, html

class LatestEntriesFeed(Feed):
    title = "Biostars.org latest"
    link = "/"
    description = "Latest 25 posts from the Biostar server"

    def items(self):
        return models.Post.objects.filter(post_type=const.POST_QUESTION).order_by('-creation_date')[:25]

    def item_title(self, item):
        return item.title

    def item_description(self, item):
        return item.html[:1000]

class LatestNewsFeed(Feed):
    title = "Biostars.org latest"
    link = "/"
    description = "Latest 25 notification for a given user"

    def title(self, obj):
        return "Entries for user %s" % obj.profile.display_name

    def get_object(self, request, uuid):
        obj = get_object_or_404(models.User, profile__uuid=uuid)
        return obj
        
    def items(self, obj):
        return models.Post.objects.filter(post_type=const.POST_QUESTION, author=obj).order_by('-creation_date')[:25]

    def item_title(self, item):
        return item.title

    def item_description(self, item):
        return item.html[:1000]