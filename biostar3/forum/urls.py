from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum.user_views import MeView
from biostar3.forum.post_views import PostList, UserList, SearchResults, PostView
from biostar3.forum.form_views import EditNode, NewAnswer, NewComment, NewPost, EditPost
from biostar3.forum import ajax

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
    url(r'^new/post/$', NewPost.as_view(), name="new_post"),
    url(r'^new/answer/(?P<pk>\d+)/$', NewAnswer.as_view(), name="new_answer"),
    url(r'^new/comment/(?P<pk>\d+)/$', NewComment.as_view(), name="new_comment"),

    # Edit existing posts.
    url(r'^edit/node/(?P<pk>\d+)/$', EditNode.as_view(), name="edit_node"),
    url(r'^edit/post/(?P<pk>\d+)/$', EditPost.as_view(), name="edit_post"),

    # Vote submission handler.
    url(r'^x/vote/$', ajax.vote_handler, name="vote_submit"),

    # template loader submission handler.
    url(r'^x/load/(?P<name>\w+)/(?P<pk>\d+)/$', ajax.load_html, name="load_html"),

)
