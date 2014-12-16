from django.shortcuts import render
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.contrib.auth import get_user_model
from . import models, query
from django.conf import settings

# Get custom user model.
User = get_user_model()


class UserList(ListView):
    """
    Generate user listing.
    """
    model = User
    template_name = "user-list.html"
    context_object_name = "users"
    paginate_by = 60

class RecentContext(ListView):
    def get_context_data(self, **kwargs):
        context = super(RecentContext, self).get_context_data(**kwargs)
        context['recent_votes'] = query.recent_votes()
        return context

class PostList(RecentContext):
    """
    Generate user listing.
    """
    model = models.Post
    template_name = "post-list.html"
    context_object_name = "posts"
    paginate_by = 60
