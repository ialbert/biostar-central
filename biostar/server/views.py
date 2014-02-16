from django.views.generic import DetailView, ListView, TemplateView
from django.conf import settings
from haystack.views import SearchView

from biostar.apps.users import auth
from biostar.apps.users.views import EditUser
from biostar.apps.posts.views import EditPost, NewPost, NewAnswer
import random

from biostar.apps.messages.models import Message
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post, Vote, Tag, Subscription
from biostar.apps.posts.auth import post_permissions
from django.contrib import messages
from django.utils.timezone import utc
from datetime import datetime, timedelta
from ordereddict import OrderedDict
from biostar import const


class BaseListMixin(ListView):
    "Base class for each mixin"
    page_title = "Title"
    paginate_by = settings.PAGINATE_BY

    def get_title(self):
        return self.page_title

    def get_context_data(self, **kwargs):
        context = super(BaseListMixin, self).get_context_data(**kwargs)
        context['page_title'] = self.get_title()

        sort =  self.request.GET.get('sort', const.POST_SORT_DEFAULT)
        limit =  self.request.GET.get('limit', const.POST_LIMIT_DEFAULT)

        if sort not in const.POST_SORT_MAP:
            messages.warning(self.request, const.POST_SORT_INVALID_MSG )
            sort = const.POST_SORT_DEFAULT

        if limit not in const.POST_LIMIT_MAP:
            messages.warning(self.request, const.POST_LIMIT_INVALID_MSG )
            limit = const.POST_LIMIT_DEFAULT

        context['sort'] = sort
        context['limit'] = limit
        return context

# The naming here needs to match that in the server_tag.py template tags.



def apply_sort(request, query):
    # Apply sort order
    sort = request.GET.get('sort', const.POST_SORT_DEFAULT)
    field = const.POST_SORT_MAP.get(sort, "-lastedit_date")
    query = query.order_by(field)

    # Apply time limit.
    limit = request.GET.get('limit', const.POST_LIMIT_DEFAULT)
    days = const.POST_LIMIT_MAP.get(limit, 0)
    if days:
        delta = const.now() - timedelta(days=days)
        query = query.filter(lastedit_date__gt=delta)
    return query


LATEST = "latest"
MYPOSTS, MYTAGS, UNANSWERED, FOLLOWING, BOOKMARKS = "myposts mytags unanswered following bookmarks".split()
POST_TYPES = dict(jobs=Post.JOB, forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)


def posts_by_topic(request, topic):
    "Returns a post query that matches a topic"
    user = request.user
    topic = topic.lower()

    if topic == MYPOSTS:
        # Get the posts that the user wrote.
        return Post.objects.my_posts(user)

    if topic == MYTAGS:
        # Get the posts that the user wrote.
        messages.success(request, 'Posts matching the <b><i class="fa fa-tag"></i> My Tags</b> setting in your user profile')
        return Post.objects.tag_search(user.profile.my_tags)

    if topic == UNANSWERED:
        # Get unanswered posts.
        return Post.objects.top_level(user).filter(type=Post.QUESTION, reply_count=0)

    if topic == FOLLOWING:
        # Get that posts that a user follows.
        messages.success(request, 'Threads that will produce notifications.')
        return Post.objects.top_level(user).filter(subs__user=user)

    if topic == BOOKMARKS:
        # Get that posts that a user bookmarked.
        return Post.objects.my_bookmarks(user)

    if topic in POST_TYPES:
        # A post type.
        return Post.objects.top_level(user).filter(type=POST_TYPES[topic])

    if topic:
        # Any type of topic.
        return Post.objects.tag_search(topic)

    # Return latest by default.
    return Post.objects.top_level(user).exclude(type=Post.BLOG)


class PostList(BaseListMixin):
    """
    This is the base class for any view that produces a list of posts.
    """
    model = Post
    template_name = "post-list.html"
    context_object_name = "posts"
    paginate_by = settings.PAGINATE_BY
    LATEST = "Latest"


    def __init__(self, *args, **kwds):
        super(PostList, self).__init__(*args, **kwds)
        self.limit = 250
        self.topic = None

    def get_title(self):
        if self.topic:
            return "%s Posts" % self.topic
        else:
            return "Latest Posts"

    def get_queryset(self):
        self.topic = self.kwargs.get("topic", "")
        query = posts_by_topic(self.request, self.topic)
        query = apply_sort(self.request, query)

        # Limit latest topics to a few pages.
        if not self.topic:
            query = query[:settings.SITE_LATEST_POST_LIMIT]
        return query

    def get_context_data(self, **kwargs):
        context = super(PostList, self).get_context_data(**kwargs)
        context['topic'] = self.topic or self.LATEST
        return context


class MessageList(ListView):
    """
    This is the base class for any view that produces a list of posts.
    """
    model = Message
    template_name = "message-list.html"
    context_object_name = "objects"
    paginate_by = settings.PAGINATE_BY

    def get_queryset(self):
        objs = Message.objects.filter(user=self.request.user).select_related("body").order_by('-creation_date')
        return objs

    def get_context_data(self, **kwargs):
        context = super(MessageList, self).get_context_data(**kwargs)
        context['topic'] = "messages"
        context['page_title'] = "Messages"
        return context


class TagList(BaseListMixin):
    """
    Produces the list of tags
    """
    model = Tag
    page_title = "Tags"
    context_object_name = "tags"
    template_name = "tag-list.html"


class VoteList(ListView):
    """
    Produces the list of votes
    """
    model = Message
    template_name = "vote-list.html"
    context_object_name = "votes"
    paginate_by = settings.PAGINATE_BY

    def get_queryset(self):
        objs = Vote.objects.filter(post__author=self.request.user).select_related("post").order_by('-date')
        return objs

    def get_context_data(self, **kwargs):
        context = super(VoteList, self).get_context_data(**kwargs)
        people = [v.author for v in context[self.context_object_name]]
        random.shuffle(people)
        context['topic'] = "votes"
        context['page_title'] = "Votes"
        context['people'] = people
        return context


class UserList(ListView):
    """
    Base class for the showing user listing.
    """
    model = User
    template_name = "user-list.html"
    context_object_name = "users"
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super(UserList, self).get_context_data(**kwargs)
        context['topic'] = "Users"
        return context


class UserDetails(DetailView):
    """
    Renders a user profile.
    """
    model = User
    template_name = "user-details.html"
    context_object_name = "target"

    def get_object(self):
        obj = super(UserDetails, self).get_object()
        obj = auth.user_permissions(request=self.request, target=obj)
        return obj

    def get_context_data(self, **kwargs):
        context = super(UserDetails, self).get_context_data(**kwargs)
        return context


class EditUser(EditUser):
    template_name = "user-edit.html"


class PostDetails(DetailView):
    """
    Shows a thread, top level post and all related content.
    """
    model = Post
    context_object_name = "post"
    template_name = "post-details.html"

    def get_object(self):
        user = self.request.user

        obj = super(PostDetails, self).get_object()

        # Adds the permissions
        obj = post_permissions(request=self.request, post=obj)

        # This will be piggybacked on the main object.
        obj.sub = Subscription.get_sub(post=obj, user=user)

        # Just a sanity check to start at top level.
        if obj != obj.root:
            obj = obj.root

        # Populate the object to build a tree that contains all posts in the thread.
        # Answers sorted before comments.
        thread = [post_permissions(request=self.request, post=post) for post in Post.objects.get_thread(obj)]

        # Do a little preprocessing.
        answers = [p for p in thread if p.type == Post.ANSWER]

        tree = OrderedDict()
        for post in thread:

            if post.type == Post.COMMENT:
                tree.setdefault(post.parent_id, []).append(post)

        store = {Vote.UP: set(), Vote.BOOKMARK: set()}

        if user.is_authenticated():
            pids = [p.id for p in thread]
            votes = Vote.objects.filter(post_id__in=pids, author=user).values_list("post_id", "type")

            for post_id, vote_type in votes:
                store.setdefault(vote_type, set()).add(post_id)

        # Shortcuts to each storage.
        bookmarks = store[Vote.BOOKMARK]
        upvotes = store[Vote.UP]

        def decorate(post):
            post.has_bookmark = post.id in bookmarks
            post.has_upvote = post.id in upvotes

        # Add attributes by mutating the objects
        map(decorate, thread + [obj])

        # Additional attributes used during rendering
        obj.tree = tree
        obj.answers = answers

        return obj

    def get_context_data(self, **kwargs):
        context = super(PostDetails, self).get_context_data(**kwargs)
        context['request'] = self.request
        return context


class SiteSearch(SearchView):
    extra_context = lambda x: dict(topic="search")


class RSS(TemplateView):
    template_name = "rss-info.html"