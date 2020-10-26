from django.contrib.syndication.views import Feed
from django.shortcuts import render

from biostar.forum.models import Post
from biostar.forum.util import now, split

from datetime import timedelta
from django.contrib.sites.models import Site

from biostar.forum.models import User,Profile
from django.conf import settings
import bleach

SITE_NAME = settings.SITE_NAME

FEED_COUNT = 25

HTML_THRESHOLD = 1500


def reduce_html(text):
    if len(text) > HTML_THRESHOLD:
        text = bleach.clean(text, strip=True)
        text = text[:HTML_THRESHOLD] + u' ... '
    return text


def info(request):
    """
    The RSS info,
    """
    admins = User.objects.filter(is_superuser=True)
    mods = User.objects.filter(profile__role=Profile.MODERATOR).exclude(id__in=admins)
    admins = admins.prefetch_related("profile").order_by("-profile__score")
    mods = mods.prefetch_related("profile").order_by("-profile__score")
    context = dict(admins=admins, mods=mods)

    return render(request, template_name="rss.html", context=context)


class PostBase(Feed):
    "Forms the base class to any feed producing posts"
    link = "/"
    title = "title"
    description = "description"

    def item_title(self, item):
        #title = item.title if item.is_question else f"{item.get_type_display()}: {item.title}"
        return item.title

    def item_description(self, item):
        return reduce_html(item.content)

    def item_guid(self, obj):
        return f"{obj.uid}"

    def item_pubdate(self, item):
        return item.creation_date


class LatestFeed(PostBase):
    "Latest posts"
    title = f"{SITE_NAME} latest!"
    description = f"Latest 25 posts from the {title}"

    def items(self):
        # Delay posts hours.
        delay_time = now() - timedelta(hours=2)
        posts = Post.objects.valid_posts(creation_date__lt=delay_time).exclude(type=Post.BLOG).order_by('-creation_date')
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
        return f"Activity on posts {text}"

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
        return f"Activity on posts {obj}"

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
        return f"Posts that match  {obj}"

    def title(self, obj):
        return "Post Feed"

    def items(self, obj):
        posts = Post.objects.filter(tags__name=obj)
        return posts[:FEED_COUNT]


class UserFeed(PostBase):

    def get_object(self, request, text):
        return text

    def description(self, obj):
        return f"Posts for users that match  {obj}"

    def title(self, obj):
        return "User Feed"

    def items(self, text):
        ids = split(text)
        posts = Post.objects.filter(author__id__in=ids).order_by('-creation_date')
        return posts[:FEED_COUNT]

