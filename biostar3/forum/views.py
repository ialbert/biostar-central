from django.shortcuts import render
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.contrib.auth import get_user_model
from .models import Post

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

class PostList(ListView):
    """
    Generate user listing.
    """
    model = Post
    template_name = "post-list.html"
    context_object_name = "posts"
    paginate_by = 60