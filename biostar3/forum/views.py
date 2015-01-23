from __future__ import absolute_import, division, print_function, unicode_literals
from django.shortcuts import render
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.contrib.auth import get_user_model
from . import models, query
from django.conf import settings
from biostar3.forum import search
from django.contrib import messages
from django import shortcuts
from django.core.urlresolvers import reverse

# Get custom user model.
User = get_user_model()


def check_request(request):
    """
    Returns a dictionary with valid sort, limit and query fields.
    """


class ExtraContext(object):
    """
    Base class that adds extra context to a list view.
    """
    html_title = "Page Title"
    page_title = "Welcome to Biostar!"

    def get_context_data(self, **kwargs):
        context = super(ExtraContext, self).get_context_data(**kwargs)
        context['recent_votes'] = query.recent_votes()
        context['html_title'] = self.html_title
        context['page_title'] = self.page_title

        # Fix potentially invalid sort and limit parameters
        sort = self.request.GET.get('sort', '')
        limit = self.request.GET.get('limit', '')

        if sort and sort not in settings.POST_SORT_MAP:
            messages.warning(self.request, settings.POST_SORT_INVALID_MSG)
            sort = ''

        context['sort'] = sort
        context['limit'] = limit
        context['q'] = self.q

        return context

    @property
    def q(self):
        """
        Search query
        """
        return self.request.GET.get('q', '')[:150]


class UserList(ExtraContext, ListView):
    """
    Generates user listing.
    """
    html_title = "Users"
    model = User
    template_name = "user_list.html"
    context_object_name = "users"
    paginate_by = 20


class PostList(ExtraContext, ListView):
    """
    Generates post lists from a web request.
    Handles filtering by tags, users and post types.
    """
    model = models.Post
    template_name = "post_list.html"
    context_object_name = "posts"
    paginate_by = settings.POST_PAGINATE_BY
    html_title = "Posts"


class SearchResults(PostList):
    """
    Handles search requests.
    """
    template_name = "post_search_results.html"
    html_title = "Search Results"

    def get_queryset(self, **kwargs):
        posts = search.plain(self.q)
        return posts

    def dispatch(self, request, *args, **kwargs):
        # check if there is some video onsite
        if not self.q:
            return shortcuts.redirect(reverse("home"))
        else:
            return super(SearchResults, self).dispatch(request, *args, **kwargs)


class PostView(ExtraContext, ListView):
    """
    Generates apage that contains a full thread.
    """
    model = models.Post
    template_name = "post_detail.html"
    context_object_name = "post"
    html_title = "Posts"
