from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum import form_views, user_views, charm_views, post_views, api, feeds, ajax

urlpatterns = patterns('',

    # List of posts.
    url(r'^$', post_views.post_list, name='home'),

    # List of posts.
    url(r'^p/list/$', post_views.post_list, name='post_list'),

    # List of unanswered questions.
    url(r'^p/unanswered/$', post_views.unanswered, name='unanswered'),

    # Renders search results.
    url(r'^site/search/$', post_views.search_results, name='search'),

    # A shortcut to a user's account.
    url(r'^site/me/$', user_views.me_view, name='me'),

    # A shortcut to a user's profile edit
    url(r'^site/edit/my/profile/$', user_views.edit_my_profile, name='edit_my_profile'),

    # The signin/signup view.
    url(r'^site/sign_up/$', user_views.sign_up, name='sign_up'),

    # This is to retarget the url from django allauth.
    url(r'^accounts/signup/$', user_views.sign_up, name='account_signup'),

    # This is to retarget the url from django allauth.
    url(r'^accounts/login/$', user_views.Login.as_view(), name='account_login'),

    # Post (thread) view.
    url(r'^p/(?P<pk>\d+)/$', post_views.post_view, name="post_view"),

    # User view.
    url(r'^u/(?P<pk>\d+)/$', user_views.user_view, name="user_view"),

    # The list of users.
    url(r'^u/list/$', user_views.user_list, name="user_list"),

    # The list of users.
    url(r'^u/edit/(?P<pk>\d+)/$', form_views.user_edit, name="user_edit"),

    # The list of badges.
    url(r'^b/list/$', user_views.badge_list, name="badge_list"),

    # The list of badges.
    url(r'^b/(?P<pk>\d+)/$', user_views.badge_view, name="badge_view"),

    # The list of users.
    url(r'^b/awards/for/(?P<pk>\d+)/$', user_views.award_list, name="award_list"),

    #
    # Site specific urls. The may be blocked via robots.txt
    #
    # Create new content: answer, comments
    url(r'^site/new/post/$', form_views.create_toplevel_post, name="new_post"),
    url(r'^site/new/page/$', form_views.create_page_post, name="new_page"),
    url(r'^site/new/answer/(?P<pk>\d+)/$', form_views.create_answer, name="new_answer"),
    url(r'^site/new/comment/(?P<pk>\d+)/$', form_views.create_comment, name="new_comment"),

    # Edit existing posts.
    url(r'^site/edit/post/(?P<pk>\d+)/$', form_views.post_edit, name="post_edit"),

    # Vote submission handler.
    url(r'^site/x/vote/$', ajax.vote_handler, name="vote_submit"),

    # Loads the comment template via ajax.
    url(r'^site/x/add_comment/(?P<pk>\d+)/$', ajax.add_comment, name="add_comment"),

    # Loads the post moderation template via ajax.
    url(r'^site/x/post_moderate/(?P<pk>\d+)/$', ajax.post_moderate, name="post_moderate"),

    # Loads the user moderation template via ajax.
    url(r'^site/x/user_moderate/(?P<pk>\d+)/$', ajax.user_moderate, name="user_moderate"),

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

    # Email handlers.
    url(r'^site/incoming/email/', ajax.email_handler, name="email_handler"),

    # Biostar charms.
    url(r'^site/charms/$', charm_views.charms, name="charm_view"),

    # Biostar charm function calls.
    url(r'^site/charms/rpc/$', charm_views.charms_rpc, name="charm_rpc"),

    #
    # Group related handlers
    #
    # Group list.
    url(r'^g/list/$', post_views.group_list, name="group_list"),

    # Group information.
    url(r'^g/info/(?P<pk>\d+)/$', post_views.group_info, name="group_info"),

    # Group login redirect. Used when redirecting subdomain login.
    url(r'^g/login/(?P<pk>\d+)/$', post_views.group_login, name="group_login"),

    # Group redirect.
    url(r'^g/redirect/(?P<pk>\d+)/$', post_views.group_redirect, name="group_redirect"),

    # Group subscriptions
    url(r'^g/subscribe/(?P<pk>\d+)/$', form_views.group_subscribe, name="group_subscribe"),

    # Group create.
    url(r'^g/create/$', form_views.group_create, name="group_create"),

    # Group edit.
    url(r'^g/edit/(?P<pk>\d+)/$', form_views.group_edit, name="group_edit"),

    # Group management.
    url(r'^g/manage/(?P<pk>\d+)/$', form_views.group_manage, name="group_manage"),

    # Group premission change.
    url(r'^g/permission/(?P<pk>\d+)/$', form_views.group_permission, name="group_permission"),

    #
    # Tag specific handlers
    #
    # Tag list.
    url(r'^t/$', post_views.tag_list, name="tag_list"),

    # Planet link.
    url(r'^v/planet/$', post_views.planet_list, name="planet_list"),


    # Filter posts by tag.
    url(r'^t/(?P<name>.+)/$', post_views.tag_filter, name="tag_filter"),

    #
    # Flatpage viewers
    #
    url(r'^page/(?P<slug>[-\w]+)/$', post_views.flatpage_view, name='page_view_main'),
    url(r'^page/(?P<domain>\w+)/(?P<slug>[-\w]+)/$', post_views.flatpage_view, name='page_view'),

    #
    # RSS feed handlers.
    #
    url(r'^feeds/latest/$', feeds.LatestFeed(), name='latest_feed'),
    url(r'^feeds/planet/$', feeds.PlanetFeed(), name='planet_feed'),
    url(r'^feeds/tag/(?P<text>\S+)/$', feeds.TagFeed(), name='tag_feed'),
    url(r'^feeds/user/(?P<text>\S+)/$', feeds.UserFeed(), name='user_feed'),
    url(r'^feeds/post/(?P<text>\S+)/$', feeds.PostFeed(), name='post_feed' ),
    url(r'^feeds/type/(?P<text>\S+)/$', feeds.PostTypeFeed(), name='post_type_feed'),

    # Api related handler.
    url(r'^api/traffic/$', api.traffic, name='api_traffic'),
    url(r'^api/user/(?P<id>\d+)/$', api.user_details, name='api_user'),
    url(r'^api/post/(?P<id>[-\w]+)/$', api.post_details, name='api_post'),
    url(r'^api/vote/(?P<id>\d+)/$', api.vote_details, name='api_vote'),
    url(r'^api/stats/day/(?P<day>\d+)/$', api.daily_stats_on_day, name='api_stats_on_day'),
    url(r'^api/stats/date/(?P<year>\d{4})/(?P<month>\d{2})/(?P<day>\d{2})/$',
        api.daily_stats_on_date, name='api_stats_on_date'),

)
