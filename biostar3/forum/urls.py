from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum.post_views import PostList, UserList, SearchResults, PostView

from biostar3.forum import form_views, user_views

from biostar3.forum import ajax

urlpatterns = patterns('',

    # List of posts.
    url(r'^$', PostList.as_view(), name='home'),

    # Renders search results.
    url(r'^site/search/$', SearchResults.as_view(), name='search'),

    # A shortcut to a user's account.
    url(r'^site/me/$', user_views.me_view, name='me'),

    # The signin/signup view.
    url(r'^site/sign_up/$', user_views.sign_up, name='sign_up'),

    # This is to disable the url from django allauth.
    url(r'^accounts/signup/$', user_views.sign_up, name='account_signup'),

    # Post details.
    url(r'^p/(?P<pk>\d+)/$', PostView.as_view(), name="post_view"),

    # The list of users.
    url(r'^user/list/$', UserList.as_view(), name="user_list"),


    # Create new content: answer, comments
    url(r'^new/post/$', form_views.create_toplevel_post, name="new_post"),
    url(r'^new/answer/(?P<parent_id>\d+)/$', form_views.create_answer, name="new_answer"),
    url(r'^new/comment/(?P<parent_id>\d+)/$', form_views.create_comment, name="new_comment"),

    # Edit existing posts.
    url(r'^edit/post/(?P<post_id>\d+)/$', form_views.edit_post, name="edit_post"),

    # Vote submission handler.
    url(r'^x/vote/$', ajax.vote_handler, name="vote_submit"),

    # Loads a template via ajax.
    url(r'^x/load/(?P<name>\w+)/(?P<pk>\d+)/$', ajax.load_html, name="load_html"),


)
