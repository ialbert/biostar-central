from __future__ import absolute_import, division, print_function, unicode_literals
import sys
from django.conf import settings
from biostar3.forum import models
from biostar3.forum.models import Post, Vote, UserGroup
from django.contrib.auth import get_user_model
from django.contrib import messages
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.utils.timezone import utc
from datetime import datetime, timedelta

User = get_user_model()

DEFAULT_GROUP = UserGroup.objects.filter(name=settings.DEFAULT_GROUP_NAME).first()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def ago(hours=0, minutes=0, days=0):
    since = now() - timedelta(days=days, hours=hours, minutes=minutes)
    return since


def positive_integer(text, upper=sys.maxint):
    try:
        value = int(text)
        value = 0 if (value < 0 or value > upper) else value
        return value

    except ValueError, exc:
        return 0


class DropDown(object):
    """
    Represents a dropdown menu with a current value and label.
    """
    choices = []  # Pairs of value/display.
    lookup = {}  # A dictionary to look up a display for a value
    default = ''  # The default value for the widget.
    order = {}  # The mapping to the order_by clause

    def __init__(self, request, value):
        # Current selection of the drop down.
        self.value = value if (value in self.lookup) else self.default
        # The label displayed in the interface.
        self.label = self.lookup.get(self.value, '???')
        if self.value != self.default:
            messages.info(request, "Sorting by: %s" % self.label)


class PostSortValidator(DropDown):
    choices = settings.POST_SORT_CHOICES
    lookup = settings.POST_SORT_MAP
    default = settings.POST_SORT_DEFAULT
    order = settings.POST_SORT_ORDER


class UserSortValidator(DropDown):
    choices = settings.USER_SORT_CHOICES
    lookup = settings.USER_SORT_MAP
    default = settings.USER_SORT_DEFAULT
    order = settings.USER_SORT_ORDER


class TimeLimitValidator(DropDown):
    choices = settings.TIME_LIMIT_CHOICES
    lookup = settings.TIME_LIMIT_MAP
    default = settings.TIME_LIMIT_DEFAULT

    def __init__(self, request, value):
        self.value = positive_integer(value, upper=10000)
        self.label = self.lookup.get(self.value, "%s days" % self.value)
        if self.value != self.default:
            messages.info(request, "Limiting to: %s" % self.label)


class PostPaginator(Paginator):
    def __init__(self, request, *args, **kwds):
        self.request = request
        self.page_num = request.GET.get('page', '1')

        sort = request.GET.get('sort', '')
        days = request.GET.get('days', '')
        self.q = request.GET.get('q', '')

        # Add the dropdowns
        self.sort = PostSortValidator(request, value=sort)
        self.days = TimeLimitValidator(request, value=days)

        super(PostPaginator, self).__init__(*args, **kwds)

    def curr_page(self):

        order_by = self.sort.order.get(self.sort.value, '?')

        # Apply the order to the object list.
        self.object_list = self.object_list.order_by(order_by)

        # Apply the time limit to the object list
        if self.days.value != self.days.default:
            # The field to look up depends on the object type
            # This should be refactored and made uniform.
            since = ago(days=self.days.value)
            if isinstance(self.sort, UserSortValidator):
                self.object_list = self.object_list.filter(profile__date_joined__gt=since)
            elif isinstance(self.sort, PostSortValidator):
                self.object_list = self.object_list.filter(creation_date__gt=since)

        try:
            pa = self.page(self.page_num)
        except PageNotAnInteger:
            # If page is not an integer, deliver first page.
            pa = self.page(1)
        except EmptyPage:
            # If page is out of range (e.g. 9999), deliver last page of results.
            pa = self.page(self.num_pages)

        # Add extra attributes to page.
        pa.q = self.q
        pa.sort = self.sort
        pa.days = self.days

        return pa


class UserPaginator(PostPaginator):
    def __init__(self, request, *args, **kwds):
        super(UserPaginator, self).__init__(request, *args, **kwds)
        sort = request.GET.get('sort', '')
        self.sort = UserSortValidator(request, sort)


def recent_votes():
    votes = Vote.objects.filter(post__status=Post.OPEN).select_related("post").order_by("-date")[
            :settings.RECENT_VOTE_COUNT]
    return votes


def get_recent_users():
    users = User.objects.all().select_related("profile").order_by("-profile__last_login")[:settings.RECENT_USER_COUNT]
    return users


def get_toplevel_posts(user, group):
    "Returns posts"
    posts = Post.objects.filter(type__in=Post.TOP_LEVEL, group=group)

    if not user.is_moderator:
        posts = posts.exclude(status=Post.DELETED)

    posts = posts.select_related("root", "author", "lastedit_user", "group").prefetch_related("tags")
    posts = posts.defer("content", "html")
    return posts


def get_all_posts(user, group):
    "Returns all posts by a user"

    posts = Post.objects.filter(author=user)

    if not user.is_moderator:
        posts = posts.exclude(status=Post.DELETED)

    posts = posts.select_related("root", "author", "lastedit_user", "group").prefetch_related("tags")
    posts = posts.defer("content", "html")

    return posts


def group_filter(self, name):
    "Performs a query to return posts that belong to a group"
    posts = Post.objects.filter(type__in=Post.TOP_LEVEL, group__name=name)
    posts = posts.select_related("root", "author", "lastedit_user").prefetch_related("tag_set")
    return posts


def get_thread(root, user):
    # Populate the object to build a tree that contains all posts in the thread.
    posts = Post.objects.filter(root=root)
    if not user.is_moderator:
        posts = posts.exclude(status=Post.DELETED)

    posts = posts.select_related("root", "author", "lastedit_user", "author__profile", "lastedit_user__profile")
    posts = posts.order_by("type", "-has_accepted", "-vote_count", "creation_date")
    posts = posts.defer("content")
    return posts