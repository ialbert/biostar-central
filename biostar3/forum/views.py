from __future__ import absolute_import, division, print_function, unicode_literals
from django.shortcuts import render
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.contrib.auth import get_user_model
from . import models, query
from django.conf import settings
from biostar3.context import modify_context

# Get custom user model.
User = get_user_model()


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

        # Add sort order, date limits and query content
        modify_context(context, self.request)

        return context


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

    def get_queryset(self, **kwargs):
        posts = query.get_posts(self.request.user)
        return posts
