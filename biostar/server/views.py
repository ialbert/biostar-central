from django.views.generic import DetailView, ListView
from django.conf import settings
from haystack.views import SearchView

from biostar.apps.users import auth
from biostar.apps.users.views import EditUser
from biostar.apps.posts.views import EditPost, NewPost, NewAnswer

from biostar.apps.users.models import User
from biostar.apps.posts.models import Post
from collections import defaultdict, OrderedDict
from biostar.apps.posts.auth import post_permissions


class PostList(ListView):
    """
    This is the base class for any view that produces a list of posts.
    """
    model = Post
    template_name = "post-list.html"
    context_object_name = "posts"
    paginate_by = settings.PAGINATE_BY
    LATEST = "Latest"
    POST_TYPE_TOPICS = dict(jobs=Post.JOB, forum=Post.FORUM, planet=Post.BLOG, pages=Post.PAGE)

    def __init__(self, *args, **kwds):
        super(PostList, self).__init__(*args, **kwds)
        self.limit = 250
        self.topic = None

    def page_title(self):
        if self.topic:
            return "%s Posts" % self.topic
        else:
            return "Latest Posts"

    def get_queryset(self):
        self.topic = self.kwargs.get("topic", "")

        # Internally topics are case insensitive.
        topic = self.topic.lower()
        if topic:
            if topic == "myposts":
                objs = Post.objects.filter(author=self.request.user)
            elif topic in self.POST_TYPE_TOPICS:
                objs = Post.objects.top_level(self.request.user).filter(type=self.POST_TYPE_TOPICS[topic])
            else:
                objs = Post.objects.top_level(self.request.user).filter(tags__name=topic).exclude(type=Post.BLOG)
        else:
            # Limit the latest posts so that engines don't crawl outside of the topics categories.
            objs = Post.objects.top_level(self.request.user).exclude(type=Post.BLOG)[:settings.SITE_LATEST_POST_LIMIT]

        return objs

    def get_context_data(self, **kwargs):
        context = super(PostList, self).get_context_data(**kwargs)
        context['topic'] = self.topic or self.LATEST
        context['page_title'] = self.page_title()
        context['show_create_button'] = True
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
        obj = super(PostDetails, self).get_object()

        # Adds the permissions
        obj = post_permissions(request=self.request, post=obj)

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

        # Add this to the object.
        obj.tree = tree
        obj.answers = answers

        return obj

    def get_context_data(self, **kwargs):
        context = super(PostDetails, self).get_context_data(**kwargs)
        context['request'] = self.request
        return context


class TopicDetails(DetailView):
    template_name = "topic-details.html"


class SiteSearch(SearchView):
    extra_context = lambda x: dict(topic="search")