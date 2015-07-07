from __future__ import unicode_literals, absolute_import, print_function

from django.contrib.syndication.views import Feed
from django.shortcuts import get_object_or_404
from biostar.apps.posts.models import Post
from biostar.apps.users.models import User
from biostar.apps.messages.models import Message
from biostar.apps.planet.models import BlogPost

from django.conf import settings
from django.contrib.sites.models import Site
from datetime import datetime, timedelta
import bleach
from biostar import const

SITE = Site.objects.get(id=settings.SITE_ID)
SITE_NAME = settings.SITE_NAME

FEED_COUNT = 25


def reduce_html(text):
    if len(text) > 1500:
        text = bleach.clean(text, strip=True)
        text = text[:1500] + u' ... '
    return text


def split(text):
    text = ''.join(text.split())
    rows = text.split('+')
    return rows

class PlanetFeed(Feed):
    "Latest posts"
    link = "/"
    FEED_COUNT = 50
    title = "%s Planet!" % SITE_NAME
    description = "Latest 50 posts of the %s" % title

    def item_title(self, item):
        try:
            title = u"%s (%s)" % (item.title, item.blog.title)
        except Exception, exc:
            title = item.title
        return title

    def item_description(self, item):
        return item.content[:250]

    def item_guid(self, obj):
        return "%s" % obj.id

    def items(self):
        posts = BlogPost.objects.select_related("blog").order_by('-creation_date')
        return posts[:FEED_COUNT]


class PostBase(Feed):
    "Forms the base class to any feed producing posts"
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
        return "%s" % obj.id

    def item_pubdate(self, item):
        return item.creation_date

class LatestFeed(PostBase):
    "Latest posts"
    title = "%s latest!" % SITE_NAME
    description = "Latest 25 posts from the %s" % title

    def items(self):
        # Delay posts hours.
        delay_time = const.now() - timedelta(hours=2)
        posts = Post.objects.filter(type__in=Post.TOP_LEVEL, status=Post.OPEN, creation_date__lt=delay_time).exclude(type=Post.BLOG).order_by('-creation_date')
        return posts[:FEED_COUNT]

class PostTypeFeed(PostBase):
    TYPE_MAP = {
        'job': Post.JOB, 'blog': Post.BLOG, 'question': Post.QUESTION,
        'forum': Post.FORUM, 'page': Post.PAGE
    }

    def get_object(self, request, text):
        words = split(text)
        codes = [self.TYPE_MAP[word] for word in words if word in self.TYPE_MAP]
        return codes, text

    def description(self, obj):
        code, text = obj
        return "Activity on posts  %s" % text

    def title(self, obj):
        return "Post Activity"

    def items(self, obj):
        codes, text = obj
        posts = Post.objects.filter(type__in=codes).order_by('-creation_date')
        return posts[:FEED_COUNT]


class PostFeed(PostBase):
    def get_object(self, request, text):
        return text

    def description(self, obj):
        return "Activity on posts  %s" % obj

    def title(self, obj):
        return "Post Activity"

    def items(self, text):
        ids = split(text)
        posts = Post.objects.filter(root_id__in=ids).order_by('-creation_date')
        return posts[:FEED_COUNT]


class TagFeed(PostBase):
    "Posts matching one or more tags"

    def get_object(self, request, text):
        elems = split(text)
        return ",".join(elems)

    def description(self, obj):
        return "Posts that match  %s" % obj

    def title(self, obj):
        return "Post Feed"

    def items(self, obj):
        posts = Post.objects.tag_search(obj)
        return posts[:FEED_COUNT]


class UserFeed(PostBase):

    def get_object(self, request, text):
        return text

    def description(self, obj):
        return "Posts for users that match  %s" % obj

    def title(self, obj):
        return "User Feed"

    def items(self, text):
        ids = split(text)
        posts = Post.objects.filter(author__id__in=ids).order_by('-creation_date')
        return posts[:FEED_COUNT]




