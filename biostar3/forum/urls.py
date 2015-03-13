from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum import form_views, user_views, post_views

from biostar3.forum import ajax

urlpatterns = patterns('',

    # List of posts.
    url(r'^$', post_views.post_list, name='home'),

    # List of posts.
    url(r'^p/list/$', post_views.post_list, name='post_list'),

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

    # The list of users.
    url(r'^u/list/$', user_views.user_list, name="user_list"),

    #
    # Site specific urls. The may be blocked via robots.txt
    #
    # Create new content: answer, comments
    url(r'^site/new/post/$', form_views.create_toplevel_post, name="new_post"),
    url(r'^site/new/answer/(?P<pk>\d+)/$', form_views.create_answer, name="new_answer"),
    url(r'^site/new/comment/(?P<pk>\d+)/$', form_views.create_comment, name="new_comment"),

    # Edit existing posts.
    url(r'^site/edit/post/(?P<pk>\d+)/$', form_views.post_edit, name="post_edit"),

    # Vote submission handler.
    url(r'^site/x/vote/$', ajax.vote_handler, name="vote_submit"),

    # Loads a template via ajax.
    url(r'^site/x/load/(?P<name>\w+)/(?P<pk>\d+)/$', ajax.load_html, name="load_html"),

    # Posts created by a user.
    url(r'^site/all/posts/created/by/u/(?P<pk>\d+)/$', post_views.posts_by_user, name="posts_by_user"),

    # Upvoted posts created by a user
    url(r'^site/upvoted/posts/created/by/u/(?P<pk>\d+)/$', post_views.upvoted_posts, name="upvoted_posts"),

    # A list of votes for a user.
    url(r'^site/all/votes/for/u/(?P<pk>\d+)/$', post_views.vote_list, name="vote_list"),

    # Posts created by a user.
    url(r'^site/my/bookmarks/$', post_views.my_bookmarks, name="my_bookmarks"),

    # User messages.
    url(r'^site/my/messages/$', post_views.my_messages, name="my_messages"),

    #
    # Group related handlers
    #
    # Group list.
    url(r'^g/list/$', post_views.group_list, name="group_list"),

    # Group redirect.
    url(r'^g/redirect/(?P<pk>\S+)/$', post_views.group_redirect, name="group_redirect"),

    # Group subscriptions
    url(r'^g/subscribe/(?P<pk>\S+)/$', form_views.group_subscribe, name="group_subscribe"),

    # Group login redirect. Used when redirecting subdomain login.
    url(r'^g/login/(?P<pk>\S+)$', post_views.group_login, name="group_login"),

    # Group create.
    url(r'^g/create/$', form_views.group_create, name="group_create"),

    # Group edit.
    url(r'^g/edit/(?P<pk>\d+)/$', form_views.group_edit, name="group_edit"),

    # Group manage.
    url(r'^g/manage/(?P<pk>\d+)/$', form_views.group_manage, name="group_manage"),

    # Group change premissions
    url(r'^g/permission/(?P<pk>\d+)/$', form_views.group_permission, name="group_permission"),


    #
    # Tag specific handlers
    #
    # Tag list.
    url(r'^t/$', post_views.tag_list, name="tag_list"),

    # Filter posts by tag.
    url(r'^t/(?P<name>\S+)/$', post_views.tag_filter, name="tag_filter"),

)
