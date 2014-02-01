from django.shortcuts import render_to_response
from django.views import generic
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from django.conf import settings


def get_int(text):
    try:
        return int(text)
    except ValueError, exc:
        return 1


def paginate(paginator, page):
    try:
        objs = paginator.page(page)
    except PageNotAnInteger:
        # If page is not an integer, deliver first page.
        objs = paginator.page(1)
    except EmptyPage:
        # If page is out of range (e.g. 9999), deliver last page of results.
        objs = paginator.page(paginator.num_pages)
    return objs


class PageBase(generic.TemplateView):
    """
    All pages will inherit the attributes of this class.
    """
    page_title = "Bioinformatics Answers on Biostars"

    def get_context_data(self, **kwargs):
        context = super(PageBase, self).get_context_data(**kwargs)
        context['page_title'] = self.page_title
        context['user'] = self.request.user

        # Populate attributes required for all views.
        self.page = context['page'] = self.request.GET.get("page", "1")
        self.query = context['query'] = self.request.GET.get("query", "")

        return context


class IndexView(PageBase):
    page_title = "Welcome to Biostars"
    template_name = "index.html"

    def get_context_data(self, **kwargs):
        LIMIT = None
        context = super(IndexView, self).get_context_data(**kwargs)
        objs = Post.objects.top_level(self.request.user)[:LIMIT]
        objs = Paginator(objs, settings.POSTS_PER_PAGE)
        posts = paginate(objs, self.page)
        context['posts'] = posts
        return context

class UserView(PageBase):
    page_title = "The Biostars"
    template_name = "users.html"