from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum.user_views import MeView
from biostar3.forum.post_views import PostList, UserList, SearchResults, PostView
from biostar3.forum.form_views import EditNode, NewNode, NewPost

urlpatterns = patterns('',

    # List of posts.
    url(r'^$', PostList.as_view(), name='home'),

    # Renders search results.
    url(r'^search/$', SearchResults.as_view(), name='search'),

    # A shortcut to a user's account.
    url(r'^me/$',MeView.as_view(), name='me'),

    # Post details.
    url(r'^p/(?P<pk>\d+)/$', PostView.as_view(), name="post_view"),

    # The list of users.
    url(r'^user/list/$', UserList.as_view(), name="user_list"),

    # Create new content: answer, comments
    url(r'^new/node/(?P<pk>\d+)/$', NewNode.as_view(), name="new_node"),
    url(r'^edit/node/(?P<pk>\d+)/$', EditNode.as_view(), name="edit_node"),

    url(r'^new/post/$', NewPost.as_view(), name="new_post"),


)
