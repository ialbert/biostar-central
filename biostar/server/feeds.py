from django.contrib.syndication.views import Feed
from django.shortcuts import get_object_or_404
from biostar.apps.posts.models import Post
from biostar.apps.users.models import User

from django.conf import settings
from django.contrib.sites.models import Site
from datetime import datetime, timedelta
import bleach

SITE = Site.objects.get(id=settings.SITE_ID)
SITE_NAME = settings.SITE_NAME

FEED_COUNT = 25


def reduce_html(text):
    if len(text) > 1500:
        text = bleach.clean(text, strip=True)
        text = text[:1500] + u' ... '
    return text


class PostBase(Feed):
    link = "/"
    title = "title"
    description = "description"

    def item_title(self, item):
        if item.type != Post.QUESTION:
            return "%s: %s" % (item.get_type_display(), item.title)
        else:
            return item.title

    def item_description(self, item):
        return reduce_html(item.content)

    def item_guid(self, obj):
        return "http://%s%s" % (SITE.domain, obj.get_absolute_url())


class LatestEntriesFeed(PostBase):
    title = "%s latest!" % SITE_NAME
    description = "Latest 25 posts from the %s" % title

    def items(self):
        posts = Post.objects.filter(type__in=Post.TOP_LEVEL).exclude(type=Post.BLOG).order_by('-creation_date')
        return posts[:FEED_COUNT]


class PostTypeFeed(PostBase):
    title = "title"

    def title(self, obj):
        return "Feed to %s" % obj

    def items(self, obj):
        obj = obj or ''
        posts = Post.objects.tag_search(obj)
        return posts[:FEED_COUNT]