from __future__ import unicode_literals, absolute_import, print_function

from django.contrib.syndication.views import Feed
from django.shortcuts import get_object_or_404
from .models import Post,  User, BlogPost, Message
from . import auth
from django.conf import settings
from django.contrib.sites.models import Site
from datetime import datetime, timedelta
import bleach


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
    """
    Latest posts from the planet"
    """
    site = Site.objects.get(id=settings.SITE_ID)
    link = "/"
    FEED_COUNT = 50
    title = "%s Planet!" % site.name
    description = "Latest 50 posts of the %s" % title

    def item_title(self, item):
        try:
            title = u"%s (%s)" % (item.title, item.blog.title)
        except Exception as exc:
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
    title = "%s latest!" % Site.objects.get_current().name
    description = "Latest 25 posts from the %s" % title

    def items(self):
        posts = Post.objects.filter(type__in=Post.TOP_LEVEL,
                                   ).exclude(type=Post.BLOG).order_by('-creation_date')
        return posts[:FEED_COUNT]


class PostTypeFeed(PostBase):
    TYPE_MAP = {
        'job': Post.JOB, 'blog': Post.BLOG, 'question': Post.QUESTION,
        'forum': Post.FORUM, 'page': Post.PAGE, 'answer': Post.ANSWER,
        'comment': Post.COMMENT,
    }

    def get_object(self, request, text):
        words = text.split("+")
        codes = [self.TYPE_MAP.get(word, 0) for word in words]
        print (codes)
        return codes, words

    def description(self, obj):
        code, words = obj
        return "Activity on posts of type: %s" % ", ".join(words)

    def title(self, obj):
        return "Post Activity"

    def items(self, obj):
        codes, words = obj
        posts = Post.objects.filter(type__in=codes,
                                    ).order_by('-creation_date')
        return posts[:FEED_COUNT]


class PostFeed(PostBase):
    def get_object(self, request, text):
        elems = text.split("+")
        return elems

    def description(self, obj):
        return "Activity on posts  %s" % ", ".join(obj)

    def title(self, obj):
        return "Post Activity"

    def items(self, obj):
        posts = Post.objects.filter(root_id__in=obj,
                                    ).order_by('-creation_date')
        return posts[:FEED_COUNT]


class TagFeed(PostBase):
    "Posts matching one or more tags"

    def get_object(self, request, text):
        elems = text.split("+")[:10]
        return elems

    def description(self, obj):
        return "Posts that match tags: %s" % ", ".join(obj)

    def title(self, obj):
        return "Post Feed"

    def items(self, obj):
        posts = Post.objects.filter(tags__name__in=obj)
        return posts[:FEED_COUNT]


class UserFeed(PostBase):
    def get_object(self, request, text):
        elems = text.split("+")[:10]
        elems = map(auth.safe_int, elems)
        return elems

    def description(self, obj):
        return "Posts for users with ids in: %s" % ", ".join(obj)

    def title(self, obj):
        return "User Feed"

    def items(self, obj):
        posts = Post.objects.filter(author__id__in=obj).order_by('-creation_date')
        return posts[:FEED_COUNT]




