from __future__ import absolute_import, division, print_function, unicode_literals
import sys
from django.conf import settings
from biostar3.forum import models
from biostar3.forum.models import Post, Vote, UserGroup
from django.contrib.auth import get_user_model
from django.contrib import messages
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

User = get_user_model()

DEFAULT_GROUP = UserGroup.objects.filter(name=settings.DEFAULT_GROUP_NAME).first()

def positve_integer(text, upper=sys.maxint):
    try:
        value = int(text)
        value = 0 if (value < 0 or value > upper) else value
        return value

    except ValueError, exc:
        return 0

class ExtendedPaginator(Paginator):
    def __init__(self, request, *args, **kwds):
        self.request = request
        self.curr_page = request.GET.get('page', '1')
        self.sort_val = request.GET.get('sort', '')
        self.q = request.GET.get('q', '')
        self.days_val = request.GET.get('days', '')

        # Protect against nonsensically la
        self.days_val = positve_integer(self.days_val, upper=10000)
        self.days_lab = ''

        if self.sort_val and self.sort_val not in settings.POST_SORT_MAP:
            messages.warning(self.request, settings.POST_SORT_INVALID_MSG)
            self.sort_val = ''

        # The label that the users sees for the sort.
        self.sort_lab = settings.POST_SORT_MAP.get(self.sort_val, "???")

        self.days_lab = settings.TIME_LIMIT_MAP.get(self.days_val, "%s days" % self.days_val)

        super(ExtendedPaginator, self).__init__(*args, **kwds)

    def get_page(self):

        try:
            pa = self.page(self.curr_page)
        except PageNotAnInteger:
            # If page is not an integer, deliver first page.
            pa = self.page(1)
        except EmptyPage:
            # If page is out of range (e.g. 9999), deliver last page of results.
            pa = self.page(self.num_pages)

        # Add extra attributes to paging to keep the correct context.
        pa.sort_val, pa.sort_lab = self.sort_val, self.sort_lab
        pa.q, pa.days_val, pa.days_lab = self.q, self.days_val, self.days_lab

        return pa


def get_page(request, object_list, per_page):
    paginator = ExtendedPaginator(request, object_list, per_page, orphans=False)
    return paginator.get_page()


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

    posts = posts.select_related("root", "author", "lastedit_user", "group").prefetch_related("tags").defer("content",
                                                                                                            "html")
    posts = posts.defer("content")
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