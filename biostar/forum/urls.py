
from django.conf import settings
from django.conf.urls.static import static
from django.urls import include, path  # For django versions from 2.0 and up
import debug_toolbar
from biostar.forum import views
from biostar.forum import ajax
import biostar.accounts.urls as account_patterns

urlpatterns = [

    # Main entry. Post listing.
    path('', views.latest, name='post_list'),
    path('votes/', views.myvotes, name='myvotes'),
    path('bookmarks/', views.bookmarks, name='bookmarks'),
    path('following/', views.following, name='following'),
    path('myposts/', views.myposts, name='myposts'),

    path('p/<str:uid>/', views.post_view, name='post_view'),

    path('new/post/', views.new_post, name='post_create'),
    path('new/answer/<str:uid>/', views.new_answer, name='post_answer'),
    path('new/comment/<str:uid>/', views.new_comment, name="create_comment"),

    path('t/(<slug:tag>/', views.tag_filter, name='tag_filter'),

    path('b/list/', views.badge_list, name='badge_list'),

    path('b/view/<str:uid>/', views.badge_view, name='badge_view'),

    #path(r'^sub/(?P<uid>[-\w]+)/$', views.subs_action, name='subs_action'),
    path('edit/post/<str:uid>/', views.edit_post, name='post_edit'),

    path('ajax/vote/', ajax.ajax_vote, name='vote'),
    path('ajax/test/', ajax.ajax_test, name='ajax_test'),
    path('ajax/subscribe/', ajax.ajax_subs, name='ajax_sub'),
    path('ajax/edit/', ajax.ajax_edit, name='ajax_edit'),

    path('moderate/<str:uid>/', views.post_moderate, name="post_moderate"),

    # Community urls
    path('community/', views.community_list, name='community_list'),

    # Include the accounts urls
    path('accounts/', include(account_patterns)),

]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)

    urlpatterns = [
        path(r'__debug__/', include(debug_toolbar.urls)),

    ] + urlpatterns




