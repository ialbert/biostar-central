from django.conf.urls import patterns, include, url
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

from django.views.generic import TemplateView
from biostar.server import views, ajax, search, moderate
from biostar.apps.posts.views import NewAnswer, NewPost, EditPost

urlpatterns = patterns('',

    # Post listing.
    url(r'^$', views.PostList.as_view(), name="home"),

    # Listing of all tags.
    url(r'^t/$', views.TagList.as_view(), name="tag-list"),

    # Topic listing.
    url(r'^t/(?P<topic>.+)/$', views.PostList.as_view(), name="topic-list"),

    # The list of users.
    url(r'^user/list/$', views.UserList.as_view(), name="user-list"),

    # User details.
    url(r'^u/(?P<pk>\d+)/$', views.UserDetails.as_view(), name="user-details"),

    # User details.
    url(r'^u/edit/(?P<pk>\d+)/$', views.EditUser.as_view(), name="user-edit"),

    # Post details.
    url(r'^p/(?P<pk>\d+)/$', views.PostDetails.as_view(), name="post-details"),

    # Change subscription view.
    url(r'^local/sub/(?P<pk>\d+)/(?P<type>\w+)/$', views.ChangeSub.as_view(), name="change-sub"),

    # A separate url for each post type.
    url(r'^p/new/post/$', NewPost.as_view(), name="new-post"),

    url(r'^p/new/answer/(?P<pid>\d+)/$', NewAnswer.as_view(post_type="answer"), name="new-answer"),
    url(r'^p/new/comment/(?P<pid>\d+)/$', NewAnswer.as_view(post_type="comment"), name="new-comment"),

    # Edit an existing post.
    url(r'^p/edit/(?P<pk>\d+)/$', EditPost.as_view(), name="post-edit"),

    # Message display.
    url(r'^local/messages/$', views.MessageList.as_view(), name="user-messages"),

    # Vote display.
    url(r'^local/votes/$', views.VoteList.as_view(), name="user-votes"),

    # Produces the moderator panel.
    url(r'^local/moderate/$', moderate.ModeratorPanel.as_view(), name="moderator-panel"),

    # Vote submission.
    url(r'^x/vote/$', ajax.vote_handler, name="vote-submit"),

    # Social login pages.
    (r'^accounts/', include('allauth.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),

    # Adding the search urls.
    url(r'^search/', search.Search.as_view(), name="search"),

     # Local robots.txt.
    url(r'^robots\.txt$', TemplateView.as_view(template_name="robots.txt", content_type='text/plain'), name='robots'),

)

from biostar.server.feeds import LatestFeed, TagFeed, UserFeed, PostFeed, PostTypeFeed

# Adding the RSS related urls.
urlpatterns += patterns('',

    # RSS info page.
    url(r'^info/rss/$', views.RSS.as_view(), name='rss'),

    # RSS feeds
    url(r'^feeds/latest/$', LatestFeed(), name='latest-feed'),

    url(r'^feeds/tag/(?P<text>[\w\-_\+]+)/$', TagFeed(), name='tag-feed'),
    url(r'^feeds/user/(?P<text>[\w\-_\+]+)/$', UserFeed(), name='user-feed'),
    url(r'^feeds/post/(?P<text>[\w\-_\+]+)/$', PostFeed(), name='post-feed' ),
    url(r'^feeds/type/(?P<text>[\w\-_\+]+)/$', PostTypeFeed(), name='post-type'),

)

# Adding the flatpages.
urlpatterns += patterns('django.contrib.flatpages.views',
    url(r'^info/about/$', 'flatpage', {'url': '/about/'}, name='about'),
    url(r'^info/faq/$', 'flatpage', {'url': '/faq/'}, name='faq'),
    url(r'^info/help/$', 'flatpage', {'url': '/about/'}, name='help'),
    url(r'^info/policy/$', 'flatpage', {'url': '/policy/'}, name='policy'),
)

# This is used only for the debug toolbar
if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
        url(r'^__debug__/', include(debug_toolbar.urls)),
    )