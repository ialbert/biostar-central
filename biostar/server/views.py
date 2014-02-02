from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.utils.encoding import smart_text
from django.conf import settings


def get_int(text):
    try:
        return int(text)
    except ValueError, exc:
        return 1


def get_page(paginator, page):
    try:
        objs = paginator.page(page)
    except PageNotAnInteger:
        # If page is not an integer, deliver first page.
        objs = paginator.page(1)
    except EmptyPage:
        # If page is out of range (e.g. 9999), deliver last page of results.
        objs = paginator.page(paginator.num_pages)
    return objs


class PageBase(TemplateView):
    """
    All pages will inherit the attributes of this class.
    """
    page_title = "Welcome"

    def get_context_data(self, **kwargs):
        context = super(PageBase, self).get_context_data(**kwargs)
        context['page_title'] = self.page_title
        context['user'] = self.request.user
        context['topic'] = None

        # Populate attributes required for all views.
        self.page = context['page'] = self.request.GET.get("page", "1")
        self.query = context['query'] = self.request.GET.get("query", "")

        return context

class Index(PageBase):
    page_title = "Welcome"
    template_name = "index.html"

    def paginate(self, objs):
        objs = Paginator(objs, settings.POSTS_PER_PAGE)
        objs = get_page(objs, self.page)
        return objs

    def get_context_data(self, **kwargs):
        context = super(Index, self).get_context_data(**kwargs)
        # We limit this so that search engines don't follow it forever
        objs = Post.objects.top_level(self.request.user)[:250]
        context['posts'] = self.paginate(objs)
        return context

class TopicList(Index):
    page_title = "Topics"
    template_name = "index.html"

    def get_context_data(self, topic):
        context = super(TopicList, self).get_context_data()
        objs = Post.objects.top_level(self.request.user).filter(tags__name=topic.lower())
        context['posts'] = self.paginate(objs)
        context['page_title'] = "%s Answers" % topic
        context['topic'] = topic
        return context

class UserList(PageBase):
    page_title = "User List"
    template_name = "user-list.html"

class UserDetails(DetailView):
    model = User
    page_title = "User Profile"
    template_name = "user-details.html"


class TagDetails(DetailView):
    page_title = "Topics"
    template_name = "tag-details.html"