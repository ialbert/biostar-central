from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum import form_views, user_views, post_views

from biostar3.forum import ajax

urlpatterns = patterns('',

    # List of posts.
    url(r'^$', post_views.post_list, name='home'),

    # Renders search results.
    url(r'^site/search/$', post_views.search_results, name='search'),

    # A shortcut to a user's account.
    url(r'^site/me/$', user_views.me_view, name='me'),

    # The signin/signup view.
    url(r'^site/sign_up/$', user_views.sign_up, name='sign_up'),

    # This is to retarget the url from django allauth.
    url(r'^accounts/signup/$', user_views.sign_up, name='account_signup'),

    #This is to retarget the url from django allauth.
    url(r'^accounts/login/$', user_views.Login.as_view(), name='account_login'),

    # Post (thread) view.
    url(r'^p/(?P<pk>\d+)/$', post_views.post_view, name="post_view"),

    # User view.
    url(r'^u/(?P<pk>\d+)/$', user_views.user_view, name="user_view"),

    # Posts created by a user.
    url(r'^site/all/posts/created/by/user/(?P<pk>\d+)/$', post_views.posts_by_user, name="posts_by_user"),

    # Upvoted posts created by a user
    url(r'^site/upvoted/posts/created/by/(?P<pk>\d+)/$', post_views.upvoted_posts, name="upvoted_posts"),

    # Posts created by a user.
    url(r'^site/my/bookmarks/$', post_views.my_bookmarks, name="my_bookmarks"),

    # Group list.
    url(r'^g/list/$', post_views.group_list, name="group_list"),

    # Group redirect.
    url(r'^g/redirect/(?P<domain>\S+)$', post_views.group_redirect, name="group_redirect"),

    # Group create.
    url(r'^g/create/$', form_views.group_create, name="group_create"),

    # Group edit.
    url(r'^g/edit/(?P<pk>\d+)/$', form_views.group_edit, name="group_edit"),

    # Tag list.
    url(r'^t/$', post_views.tag_list, name="tag_list"),

    # Filter posts by tag.
    url(r'^t/(?P<name>\S+)/$', post_views.tag_filter, name="tag_filter"),


    # The list of users.
    url(r'^user/list/$', user_views.user_list, name="user_list"),


    # Create new content: answer, comments
    url(r'^new/post/$', form_views.create_toplevel_post, name="new_post"),
    url(r'^new/answer/(?P<parent_id>\d+)/$', form_views.create_answer, name="new_answer"),
    url(r'^new/comment/(?P<parent_id>\d+)/$', form_views.create_comment, name="new_comment"),

    # Edit existing posts.
    url(r'^edit/post/(?P<pk>\d+)/$', form_views.edit_post, name="edit_post"),

    # Vote submission handler.
    url(r'^x/vote/$', ajax.vote_handler, name="vote_submit"),

    # Loads a template via ajax.
    url(r'^x/load/(?P<name>\w+)/(?P<pk>\d+)/$', ajax.load_html, name="load_html"),


)
