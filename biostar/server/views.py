from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView, ListView
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.utils.encoding import smart_text
from django.conf import settings
from haystack.views import SearchView


class PostList(ListView):
    model = Post
    template_name = "post-list.html"
    context_object_name = "posts"
    paginate_by = 25
    LATEST = "Latest"

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
        self.topic = self.kwargs.get("topic")
        if self.topic:
            objs = Post.objects.top_level(self.request.user).filter(tags__name=self.topic.lower())
        else:
            # Limit the latest posts so that engines don't crawl outside of the topics catergories.
            objs = Post.objects.top_level(self.request.user)[:self.limit]

        return objs

    def get_context_data(self, **kwargs):
        context = super(PostList, self).get_context_data(**kwargs)
        context['topic'] = self.topic or self.LATEST
        context['page_title'] = self.page_title()
        return context


class UserList(ListView):
    model = User
    template_name = "user-list.html"
    context_object_name = "users"
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super(UserList, self).get_context_data(**kwargs)
        context['topic'] = "Users"
        return context

class UserDetails(DetailView):
    model = User
    template_name = "user-details.html"

class PostDetails(DetailView):
    model = Post
    template_name = "post-details.html"

class TopicDetails(DetailView):
    template_name = "topic-details.html"

class SiteSearch(SearchView):
    extra_context = lambda x: dict(topic="search")