from django.views.generic import DetailView, ListView, TemplateView, RedirectView, View
from django.conf import settings

from biostar.apps.users import auth
from biostar.apps.users.views import EditUser
import random

from biostar.apps.messages.models import Message
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post, Vote, Tag, Subscription
from biostar.apps.posts.auth import post_permissions
from django.contrib import messages
from datetime import datetime, timedelta
from biostar.const import OrderedDict
from biostar import const
from braces.views import LoginRequiredMixin
from django import shortcuts
from django.http import HttpResponseRedirect
from django.core.paginator import Paginator

class BaseListMixin(ListView):
    "Base class for each mixin"
    page_title = "Title"
    paginate_by = settings.PAGINATE_BY

    def get_title(self):
        return self.page_title

    def get_context_data(self, **kwargs):
        context = super(BaseListMixin, self).get_context_data(**kwargs)
        context['page_title'] = self.get_title()

        sort = self.request.GET.get('sort', const.POST_SORT_DEFAULT)
        limit = self.request.GET.get('limit', const.POST_LIMIT_DEFAULT)

        if sort not in const.POST_SORT_MAP:
            messages.warning(self.request, const.POST_SORT_INVALID_MSG)
            sort = const.POST_SORT_DEFAULT

        if limit not in const.POST_LIMIT_MAP:
            messages.warning(self.request, const.POST_LIMIT_INVALID_MSG)
            limit = const.POST_LIMIT_DEFAULT

        context['sort'] = sort
        context['limit'] = limit
        context['q'] = self.request.GET.get('q', '')

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
POST_TYPES = dict(job=Post.JOB, forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)

def posts_by_topic(request, topic):
    "Returns a post query that matches a topic"
    user = request.user
    topic = topic.lower()

    if topic == MYPOSTS:
        # Get the posts that the user wrote.
        return Post.objects.my_posts(user)

    if topic == MYTAGS:
        # Get the posts that the user wrote.
        messages.success(request,
                         'Posts matching the <b><i class="fa fa-tag"></i> My Tags</b> setting in your user profile')
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

    if topic and topic != LATEST:
        # Any type of topic.
        return Post.objects.tag_search(topic)

    # Return latest by default.
    return Post.objects.top_level(user).exclude(type=Post.BLOG)

def reset_counts(request, label):
    "Resets counts in the session"
    label = label.lower()
    counts = request.session.get(settings.SESSION_KEY)
    if label in counts:
        counts[label] = ''
        request.session[settings.SESSION_KEY] = counts

class PostList(BaseListMixin):
    """
    This is the base class for any view that produces a list of posts.
    """
    model = Post
    template_name = "post_list.html"
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
        session = self.request.session

        context = super(PostList, self).get_context_data(**kwargs)
        context['topic'] = self.topic or self.LATEST

        reset_counts(self.request, self.topic)

        return context


class MessageList(ListView):
    """
    This is the base class for any view that produces a list of posts.
    """
    model = Message
    template_name = "message_list.html"
    context_object_name = "objects"
    paginate_by = settings.PAGINATE_BY
    topic = "messages"

    def get_queryset(self):
        objs = Message.objects.filter(user=self.request.user).select_related("body").order_by('-sent_at')
        return objs

    def get_context_data(self, **kwargs):
        context = super(MessageList, self).get_context_data(**kwargs)
        people = [m.body.author for m in context[self.context_object_name]]
        context['topic'] = self.topic
        context['page_title'] = "Messages"
        context['people'] = people
        reset_counts(self.request, self.topic)
        return context


class TagList(BaseListMixin):
    """
    Produces the list of tags
    """
    model = Tag
    page_title = "Tags"
    context_object_name = "tags"
    template_name = "tag_list.html"


class VoteList(ListView):
    """
    Produces the list of votes
    """
    model = Message
    template_name = "vote_list.html"
    context_object_name = "votes"
    paginate_by = settings.PAGINATE_BY
    topic = "votes"

    def get_queryset(self):
        objs = Vote.objects.filter(post__author=self.request.user).select_related("post").order_by('-date')
        return objs

    def get_context_data(self, **kwargs):
        context = super(VoteList, self).get_context_data(**kwargs)
        people = [v.author for v in context[self.context_object_name]]
        random.shuffle(people)
        context['topic'] = self.topic
        context['page_title'] = "Votes"
        context['people'] = people
        reset_counts(self.request, self.topic)
        return context


class UserList(ListView):
    """
    Base class for the showing user listing.
    """
    model = User
    template_name = "user_list.html"
    context_object_name = "users"
    paginate_by = 60

    def get_queryset(self):

        self.q = self.request.GET.get('q', '')
        self.sort = self.request.GET.get('sort', const.USER_SORT_DEFAULT)
        self.limit = self.request.GET.get('limit', const.POST_LIMIT_DEFAULT)

        # Apply the sort on users
        obj = User.objects.get_users(sort=self.sort, limit=self.limit, q=self.q)
        return obj

    def get_context_data(self, **kwargs):
        context = super(UserList, self).get_context_data(**kwargs)
        context['topic'] = "Users"

        context['sort'] = self.sort
        context['limit'] = self.limit
        context['q'] = self.q

        return context

class BaseDetailMixin(DetailView):

    def get_context_data(self, **kwargs):
        context = super(BaseDetailMixin, self).get_context_data(**kwargs)
        sort = self.request.GET.get('sort', const.POST_SORT_DEFAULT)
        limit = self.request.GET.get('limit', const.POST_LIMIT_DEFAULT)

        context['sort'] = sort
        context['limit'] = limit
        context['q'] = self.request.GET.get('q', '')
        return context

class UserDetails(BaseDetailMixin):
    """
    Renders a user profile.
    """
    model = User
    template_name = "user_details.html"
    context_object_name = "target"

    def get_object(self):
        obj = super(UserDetails, self).get_object()
        obj = auth.user_permissions(request=self.request, target=obj)
        return obj

    def get_context_data(self, **kwargs):
        context = super(UserDetails, self).get_context_data(**kwargs)
        target = context[self.context_object_name]
        posts = Post.objects.filter(author=target).defer("content").order_by("-creation_date")
        paginator = Paginator(posts, 10)
        page = self.request.GET.get("page", 1)
        page_obj = paginator.page(int(page))
        context['page_obj'] = page_obj
        context['posts'] = page_obj.object_list

        return context


class EditUser(EditUser):
    template_name = "user_edit.html"


class PostDetails(DetailView):
    """
    Shows a thread, top level post and all related content.
    """
    model = Post
    context_object_name = "post"
    template_name = "post_details.html"

    def get(self, *args, **kwargs):
        # This will scroll the page to the right anchor.
        self.object = self.get_object()
        context = self.get_context_data(object=self.object)

        if not self.object.is_toplevel:
            return HttpResponseRedirect(self.object.get_absolute_url())

        return self.render_to_response(context)

    def get_object(self):
        user = self.request.user

        obj = super(PostDetails, self).get_object()

        # Update the post views.
        Post.update_post_views(obj, request=self.request)

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


class ChangeSub(LoginRequiredMixin, View):
    pk, type = 0, 0
    TYPE_MAP = {"local": const.LOCAL_MESSAGE, "email": const.EMAIL_MESSAGE}

    def get(self, *args, **kwargs):
        # TODO needs to be done via POST.
        pk = self.kwargs["pk"]
        new_type = self.kwargs["type"]

        new_type = self.TYPE_MAP.get(new_type, None)

        user = self.request.user
        post = Post.objects.get(pk=pk)

        subs = Subscription.objects.filter(post=post, user=user)
        if new_type is None:
            subs.delete()
        else:
            if subs:
                subs.update(type=new_type)
            else:
                Subscription.objects.create(post=post, user=user, type=new_type)

        return shortcuts.redirect(post.get_absolute_url())


class RSS(TemplateView):
    template_name = "rss_info.html"