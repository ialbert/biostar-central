from __future__ import absolute_import, division, print_function, unicode_literals

# Python modules.
from collections import OrderedDict, defaultdict

# Django specific modules.
from django.shortcuts import render, redirect, render_to_response
from django.views.generic import DetailView, ListView, TemplateView, UpdateView, View
from django.template import RequestContext
from django.contrib.auth import get_user_model
from django.conf import settings
from django.contrib import messages
from django.core.urlresolvers import reverse

# Biostar specific local modules.
from . import models, query, search

# Get custom user model.
User = get_user_model()


class ExtraContext(object):
    """
    Base class that adds extra context to a class based view.

    Reads request URL parameters such as sort, limit and q (query) and adds them to the page
    context thus makes them available to the templates.
    """
    html_title = "Page Title"
    page_title = "Welcome to Biostar!"

    def get_context_data(self, **kwargs):

        context = super(ExtraContext, self).get_context_data(**kwargs)
        context['recent_votes'] = query.recent_votes()
        context['html_title'] = self.html_title
        context['page_title'] = self.page_title

        # Add empty values if not present in the request.
        sort = self.request.GET.get('sort', '')
        limit = self.request.GET.get('limit', '')
        q = self.request.GET.get('q', '')[:150]

        # Flash a warning message on invalid parameters.
        if sort and sort not in settings.POST_SORT_MAP:
            messages.warning(self.request, settings.POST_SORT_INVALID_MSG)
            sort = ''

        context['sort'] = sort
        context['limit'] = limit
        context['q'] = q

        return context


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

        user = self.request.user
        self.object = self.get_object()

        # This will redirect to top level and scroll the page to the right anchor.
        if not self.object.is_toplevel:
            return redirect(self.object.get_absolute_url())

        # Gets all objects in a thread. Moderators get deleted objects as well.
        thread = [p for p in query.get_thread(self.object, user)]

        # Store answers in a separate list for simpler access.
        self.object.answers = filter(lambda p: p.type == models.Post.ANSWER, thread)

        # Comments will be stored in a dictionary for fast access.
        comment_list = filter(lambda p: p.type == models.Post.COMMENT, thread)

        # Collect comments into a dictionary keyed by the parent with the posts as list.
        self.object.comments = OrderedDict()
        for post in comment_list:
            self.object.comments.setdefault(post.parent.id, []).append(post)

        # Add oject to the context.
        context = self.get_context_data(object=self.object)

        return self.render_to_response(context)


def ratelimited(request, exc):
    data = dict()
    context_instance = RequestContext(request)
    return render_to_response("ratelimit_info.html", data, context_instance)
