from __future__ import absolute_import, division, print_function, unicode_literals
from django.shortcuts import render
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.contrib.auth import get_user_model
from . import models, query
from django.conf import settings
from biostar3.forum import search
from django.contrib import messages

# Get custom user model.
User = get_user_model()


def check_request(request):
    """
    Returns a dictionary with valid sort, limit and query fields.
    """


class ExtraContext(ListView):
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


class UserList(ExtraContext):
    """
    Generate user listing.
    """
    html_title = "Users"
    model = User
    template_name = "user_list.html"
    context_object_name = "users"
    paginate_by = 60


class PostList(ExtraContext):
    """
    Generate post lists from a request.
    """
    model = models.Post
    template_name = "post_list.html"
    context_object_name = "posts"
    paginate_by = settings.POST_PAGINATE_BY
    html_title = "Posts"


class SearchResults(PostList):
    template_name = "post_search_results.html"
    html_title = "Search Results"

    def get_queryset(self, **kwargs):
        posts = search.plain(self.q)
        return posts