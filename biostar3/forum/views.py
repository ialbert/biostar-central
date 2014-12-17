from django.shortcuts import render
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.contrib.auth import get_user_model
from . import models, query
from django.conf import settings

# Get custom user model.
User = get_user_model()


class UserList(ListView):
    """Generate user listing."""
    model = User
    template_name = "user-list.html"
    context_object_name = "users"
    paginate_by = 60


class ExtraContext(ListView):
    """Base class that add extra context."""

    def get_context_data(self, **kwargs):
        context = super(ExtraContext, self).get_context_data(**kwargs)
        context['recent_votes'] = query.recent_votes()
        return context


class PostList(ExtraContext):
    """Generate post lists from a request."""
    model = models.Post
    template_name = "post-list.html"
    context_object_name = "posts"
    paginate_by = settings.POST_PAGINATE_BY

    def get_queryset(self, **kwargs):
        posts = query.get_posts(self.request.user)
        return posts