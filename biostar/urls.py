from django.conf.urls import patterns, include, url
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

from django.views.generic import TemplateView
from biostar.server import views, ajax, search, moderate, api
from biostar.apps.posts.views import NewAnswer, NewPost, EditPost, external_post_handler
from biostar.apps.users.views import external_logout, external_login, CaptchaView, EmailListView
from biostar.apps.planet.views import BlogPostList

urlpatterns = patterns('',

    # Post listing.
    url(r'^$', views.PostList.as_view(), name="home"),

    # Listing of all tags.
    url(r'^t/$', views.TagList.as_view(), name="tag-list"),

    # Badge view details.
    url(r'^b/(?P<pk>\d+)/$', views.BadgeView.as_view(), name="badge-view"),

    # Badge list details.
    url(r'^b/list/$', views.BadgeList.as_view(), name="badge-list"),

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
    url(r'^p/new/post/$', views.RateLimitedNewPost.as_view(), name="new-post"),

    # A new external post
    url(r'^p/new/external/post/$', external_post_handler, name="new-external-post"),

    url(r'^p/new/answer/(?P<pid>\d+)/$', views.RateLimitedNewAnswer.as_view(post_type="answer"), name="new-answer"),
    url(r'^p/new/comment/(?P<pid>\d+)/$', views.RateLimitedNewAnswer.as_view(post_type="comment"), name="new-comment"),

    # Edit an existing post.
    url(r'^p/edit/(?P<pk>\d+)/$', EditPost.as_view(), name="post-edit"),

    # Message display.
    url(r'^local/list/$', EmailListView.as_view(), name="email-list"),

    # Message display.
    url(r'^local/messages/$', views.MessageList.as_view(), name="user-messages"),

    # Vote display.
    url(r'^local/votes/$', views.VoteList.as_view(), name="user-votes"),

    # Produces the moderator panel.
    url(r'^local/moderate/post/(?P<pk>\d+)/$', moderate.PostModeration.as_view(), name="post-moderation"),

    # Produces the moderator panel.
    url(r'^local/moderate/user/(?P<pk>\d+)/$', moderate.UserModeration.as_view(), name="user-moderation"),

    # Full login and logout
    url(r'^site/login/$', external_login, name="login"),
    url(r'^site/logout/$', external_logout, name="logout"),
    url(r'^accounts/signup/$', CaptchaView.as_view(), name="signup"),

    # Email handlers
    url(r'^local/email/', views.email_handler, name="email-handler"),

    # Search the body.
    url(r'^local/search/page/', search.Search.as_view(), name="search-page"),

    # Search the titles.
    url(r'^local/search/title/', search.search_title, name="search-title"),

    # Returns suggested tags
    url(r'^local/search/tags/', search.suggest_tags, name="suggest-tags"),


    # Returns the planet view
    url(r'^planet/$', BlogPostList.as_view(), name="planet"),

    # Vote submission.
    url(r'^x/vote/$', ajax.vote_handler, name="vote-submit"),

    # Social login pages.
    (r'^accounts/', include('allauth.urls')),

    # Redirecting old posts urls from previous versions of Biostar
    url(r'^post/redirect/(?P<pid>\d+)/$', views.post_redirect),
    url(r'^post/show/(?P<pid>\d+)/$', views.post_redirect),
    url(r'^post/show/(?P<pid>\d+)/([-\w]+)/$', views.post_redirect),
    url(r'^questions/(?P<pid>\d+)/$', views.post_remap_redirect),
    url(r'^questions/(?P<pid>\d+)/([-\w]+)/$', views.post_remap_redirect),
    url(r'^questions/tagged/(?P<tag>.+)/$',views.tag_redirect),

    # Api.
    url(r'^api/ping/$', api.ping, name='api-ping'),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),

     # Local robots.txt.
    url(r'^robots\.txt$', TemplateView.as_view(template_name="robots.txt", content_type='text/plain'), name='robots'),

)

# Adding the sitemap.
urlpatterns += patterns('',
    (r'^sitemap\.xml$', 'django.contrib.sitemaps.views.sitemap', {'sitemaps': search.sitemaps})
)

from biostar.server.feeds import LatestFeed, TagFeed, UserFeed, PostFeed, PostTypeFeed, PlanetFeed

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
    url(r'^feeds/planet/$', PlanetFeed(), name='planet-feed'),
)

urlpatterns += patterns('',
    url(r'^info/(?P<slug>\w+)/$', views.FlatPageView.as_view(), name='flatpage'),
    url(r'^info/update/(?P<pk>\d+)/$', views.FlatPageUpdate.as_view(), name='flatpage-update'),
)

# This is used only for the debug toolbar
if settings.DEBUG:
    import debug_toolbar
    from biostar.apps.users.views import test_login
    urlpatterns += patterns('',
        url(r'^__debug__/', include(debug_toolbar.urls)),
        url(r'^test/login/', test_login),
    )