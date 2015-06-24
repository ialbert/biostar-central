from __future__ import absolute_import, division, print_function, unicode_literals
import sys
from django.conf import settings
from biostar3.forum import models
from biostar3.forum.models import Post, Vote, Award
from django.contrib.auth import get_user_model
from django.contrib import messages
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.utils.timezone import utc
from datetime import datetime, timedelta

User = get_user_model()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def ago(hours=0, minutes=0, days=0):
    since = now() - timedelta(days=days, hours=hours, minutes=minutes)
    return since


def positive_integer(text, upper=10 ^ 9):
    try:
        value = int(text)
        value = 0 if (value < 0 or value > upper) else value
        return value

    except ValueError as exc:
        return 0


class DropDown(object):
    """
    Represents a dropdown menu with a current value and label.
    """
    # Pairs of value/display.
    choices = []
    # A dictionary to look up a display for a value
    lookup = {}
    # The default value for the widget.
    default = ''
    # The mapping to the order_by clause
    order = {}

    def __init__(self, request, value):
        # Current selection of the drop down.
        self.value = value if (value in self.lookup) else self.default
        # The label displayed in the interface.
        self.label = self.lookup.get(self.value, '')
        if self.value != self.default:
            messages.info(request, "Sorting by: %s" % self.label)

    def time_filter(self, queryset, value):
        return queryset


class PostSortValidator(DropDown):
    choices = settings.POST_SORT_CHOICES
    lookup = settings.POST_SORT_MAP
    default = settings.POST_SORT_DEFAULT
    order = settings.POST_SORT_ORDER

    def time_filter(self, queryset, value):
        return queryset.filter(creation_date__gt=value)


class UserSortValidator(DropDown):
    choices = settings.USER_SORT_CHOICES
    lookup = settings.USER_SORT_MAP
    default = settings.USER_SORT_DEFAULT
    order = settings.USER_SORT_ORDER

    def time_filter(self, queryset, value):
        return queryset.filter(profile__date_joined__gt=value)


class TagSortValidator(DropDown):
    choices = [("asc", "Alphabetical"), ("desc", "Reversed")]
    lookup = dict(choices)
    default = "asc"
    order = dict(asc="name", desc="-name")


class GroupSortValidator(DropDown):
    choices = [("asc", "Alphabetical"), ("desc", "Reversed")]
    lookup = dict(choices)
    default = "asc"
    order = dict(asc="name", desc="-name")


class TimeLimitValidator(DropDown):
    choices = settings.TIME_LIMIT_CHOICES
    lookup = settings.TIME_LIMIT_MAP
    default = settings.TIME_LIMIT_DEFAULT

    def __init__(self, request, value):
        self.value = positive_integer(value, upper=10000)
        self.label = self.lookup.get(self.value, "%s days" % self.value)
        if self.value != self.default:
            messages.info(request, "Limiting to: %s" % self.label)


class ExtendedPaginator(Paginator):
    def __init__(self, request, sort_class=DropDown, time_class=DropDown, *args, **kwds):

        self.request = request
        self.page_num = request.GET.get('page', '1')

        sort = request.GET.get('sort', '')
        days = request.GET.get('days', '')
        self.q = request.GET.get('q', '')

        # Add the dropdowns
        self.sort = sort_class(request, value=sort)
        self.days = time_class(request, value=days)

        super(ExtendedPaginator, self).__init__(*args, **kwds)

    def curr_page(self):

        order_by = self.sort.order.get(self.sort.value, '')

        # Apply the order to the object list.
        if self.sort.order:
            self.object_list = self.object_list.order_by(order_by)

        # Apply the time limit to the object list
        if self.days.value != self.days.default:
            # The field to look up depends on the object type
            since = ago(days=self.days.value)
            self.object_list = self.sort.time_filter(self.object_list, since)

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


def recent_votes(request):
    votes = Vote.objects.filter(post__status=Post.OPEN) \
                .select_related("post", "post__root", "post__author").order_by("-date")[:settings.RECENT_VOTE_COUNT]
    return votes


def recent_users(request):
    users = User.objects.exclude(status=User.BANNED).select_related("profile") \
                .order_by("-profile__last_login")[:settings.RECENT_USER_COUNT]
    return users


def recent_awards(request):
    users = Award.objects.filter().select_related("post", "user") \
                .order_by("-date")[:settings.RECENT_AWARD_COUNT]
    return users


def recent_replies(request):
    posts = Post.objects.filter().select_related("author", "post", "root").exclude(type__in=Post.TOP_LEVEL) \
                .order_by("-creation_date")[:settings.RECENT_USER_COUNT]
    return posts



def get_toplevel_posts(user):
    "Returns posts"
    posts = Post.objects.filter(type__in=Post.QUERY_TYPES)

    if not user.is_moderator:
        posts = posts.exclude(status=Post.DELETED)

    posts = posts.select_related("root", "author", "lastedit_user").prefetch_related("tags")
    posts = posts.defer("content", "html")
    return posts


def get_all_posts(request, target):
    "Returns all posts by a user"

    posts = Post.objects.filter(author=target)

    if not request.user.is_moderator:
        posts = posts.exclude(status=Post.DELETED)

    posts = posts.select_related("root", "author", "lastedit_user").prefetch_related("tags")
    posts = posts.defer("content", "html")

    return posts


def get_posts_by_vote(user,  vote_types):
    posts = Post.objects.filter(votes__type__in=vote_types, votes__post__author=user)
    posts = posts.distinct()
    return posts


def get_my_bookmarks(user):
    posts = Post.objects.filter(votes__type=Vote.BOOKMARK, votes__author=user)
    posts = posts.distinct()
    return posts


def site_filter(self, sites=[]):
    "Performs a query to return posts that belong to a group"
    posts = Post.objects.filter(type__in=Post.TOP_LEVEL, site_id__in=sites)
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