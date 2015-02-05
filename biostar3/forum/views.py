from __future__ import absolute_import, division, print_function, unicode_literals
from django.shortcuts import render, redirect
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.contrib.auth import get_user_model
from . import models, query
from django.conf import settings
from biostar3.forum import search
from django.contrib import messages
from django.core.urlresolvers import reverse
from collections import OrderedDict, defaultdict


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

    def get_queryset(self):
        results = query.get_toplevel_posts(user=self.request.user, group=self.request.group)
        return results


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
            return redirect(reverse("home"))
        else:
            return super(SearchResults, self).dispatch(request, *args, **kwargs)


class PostView(ExtraContext, DetailView):
    """
    Generates apage that contains a full thread.
    """
    model = models.Post
    template_name = "post_detail.html"
    context_object_name = "post"
    html_title = "Posts"

    def get(self, *args, **kwargs):
        # This will scroll the page to the right anchor.
        self.object = self.get_object()
        context = self.get_context_data(object=self.object)

        if not self.object.is_toplevel:
            return redirect(self.object.get_absolute_url())

        return self.render_to_response(context)

    def get_object(self, *args, **kwargs):

        user = self.request.user
        root = super(PostView, self).get_object()

        # The correct data representation would be an ordered tree amp.
        thread = [p for p in query.get_thread(root, user)]

        # Store answers in a separate list.
        answers = filter(lambda p: p.type == models.Post.ANSWER, thread)

        # Populate the answers.
        root.answers = answers

        # Comments will be stored in a dictionary for fast access.
        comment_list = filter(lambda p: p.type == models.Post.COMMENT, thread)

        comments = OrderedDict()
        for post in comment_list:
            comments.setdefault(post.id, []).append(post)

        # Populate comments.
        root.comments = comments

        return root