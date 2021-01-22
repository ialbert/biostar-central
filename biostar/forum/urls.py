from django.conf import settings
from django.contrib import admin
from django.conf.urls.static import static
from django.urls import include, path  # For django versions from 2.0 and up
import debug_toolbar
from biostar.forum import views
from biostar.accounts.views import image_upload_view
from biostar.forum import ajax, api, feed
from biostar.accounts.urls import account_patterns

from biostar.planet.urls import planet_patterns


forum_patterns = [

    # Main entry. Post listing.
    path('', views.latest, name='post_list'),

    path(r'pages/<str:fname>/', views.pages, name='pages'),

    path('votes/', views.myvotes, name='myvotes'),
    path('bookmarks/', views.bookmarks, name='bookmarks'),
    path('following/', views.following, name='following'),
    path('myposts/', views.myposts, name='myposts'),
    path('mytags/', views.mytags, name='mytags'),
    path('p/<str:uid>/', views.post_view, name='post_view'),
    path('post/search/', views.post_search, name='post_search'),

    path('new/post/', views.new_post, name='post_create'),

    path('b/list/', views.badge_list, name='badge_list'),
    path('t/list/', views.tags_list, name='tags_list'),

    path('b/view/<str:uid>/', views.badge_view, name='badge_view'),

    # Ajax calls
    path('ajax/digest/', ajax.ajax_digest, name='ajax_digest'),
    path('ajax/vote/', ajax.ajax_vote, name='vote'),
    path('ajax/test/', ajax.ajax_test, name='ajax_test'),
    path('ajax/subscribe/', ajax.ajax_subs, name='ajax_sub'),
    path('ajax/delete/', ajax.ajax_delete, name='ajax_delete'),
    path('drag/and/drop/', ajax.drag_and_drop, name='drag_and_drop'),
    path('similar/posts/<str:uid>/', ajax.similar_posts, name='similar_posts'),
    path('ajax/digest/<str:uid>/', ajax.ajax_digest, name='ajax_digest'),
    path('ajax/edit/<str:uid>/', ajax.ajax_edit, name='ajax_edit'),
    path('ajax/comment/create/', ajax.ajax_comment_create, name='ajax_comment_create'),
    path('inplace/form/', ajax.inplace_form, name='inplace_form'),
    path('ajax/user/image/<str:username>/', ajax.user_image, name='user_image'),
    path('similar/posts/<str:uid>/', ajax.similar_posts, name='similar_posts'),
    path('ajax/report/spam/<str:post_uid>/', ajax.report_spam, name='report_spam'),
    path('release/<str:uid>/', ajax.release_suspect, name='release_suspect'),
    #path('most/recent/users/', ajax.most_recent_users, name='most_recent_users'),

    path('moderate/<str:uid>/', views.post_moderate, name="post_moderate"),

    path(r'mark/spam/<str:uid>/', views.mark_spam, name='mark_spam'),
    path(r'mark/spam/<str:uid>/', views.release_quar, name='release_quar'),

    # Community urls
    path('community/', views.community_list, name='community_list'),
    path('ajax/handle/search/', ajax.handle_search, name='handle_search'),

    # Api calls
    path(r'api/traffic/', api.traffic, name='api_traffic'),
    path(r'api/user/<str:uid>/', api.user_details, name='api_user'),
    path(r'api/tags/list/', api.tags_list, name='api_tags_list'),
    path(r'api/post/<str:uid>/', api.post_details, name='api_post'),
    path(r'api/vote/<str:uid>/', api.vote_details, name='api_vote'),
    path(r'api/watched/tags/<str:email>/', api.watched_tags, name='api_tags'),
    path(r'api/email/<str:email>/', api.user_email, name='user_email'),
    path(r'api/stats/day/<int:day>/', api.daily_stats_on_day, name='api_stats_on_day'),
    path(r'api/stats/date/<int:year>/<int:month>/<int:day>/', api.daily_stats_on_date,
         name='api_stats_on_date'),


    # RSS info page.
    path(r'info/rss/', feed.info, name='rss'),

    # RSS feeds
    path(r'feeds/latest/', feed.LatestFeed(), name='latest_feed'),
    path(r'feeds/tag/<str:text>/', feed.TagFeed(), name='tag_feed'),
    path(r'feeds/user/<str:text>/', feed.UserFeed(), name='user_feed'),
    path(r'feeds/post/<str:text>/', feed.PostFeed(), name='post_feed' ),
    path(r'feeds/type/<str:text>/', feed.PostTypeFeed(), name='post_type'),
    #path(r'^feeds/planet/$', feed.PlanetFeed(), name='planet-feed'),

]


urlpatterns = [

    # Main entry. Post listing.
    path('', include(forum_patterns)),

    # Include the accounts urls
    path('accounts/', include(account_patterns)),

    # Include the planet urls
    path('planet/', include(planet_patterns)),

]

if settings.PAGEDOWN_IMAGE_UPLOAD_ENABLED:

    urlpatterns += [
        # Pagedown image upload url.
        path('pagedown/image-upload/', image_upload_view, name="pagedown-image-upload")
    ]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)

if settings.DEBUG and settings.DEBUG_TOOLBAR:
    urlpatterns += [
          path('__debug__/', include(debug_toolbar.urls)),
    ]
