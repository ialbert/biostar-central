
from django.conf import settings
from django.contrib import admin
from django.conf.urls.static import static
from django.urls import include, path  # For django versions from 2.0 and up
import debug_toolbar
from biostar.forum import views
from biostar.forum import ajax
import biostar.accounts.urls as account_patterns

urlpatterns = [

    # Main entry. Post listing.
    path('', views.latest, name='post_list'),

    path('policy/', views.policy, name='policy'),

    path('votes/', views.myvotes, name='myvotes'),
    path('bookmarks/', views.bookmarks, name='bookmarks'),
    path('following/', views.following, name='following'),
    path('myposts/', views.myposts, name='myposts'),
    path('mytags/', views.mytags, name='mytags'),
    path('p/<str:uid>/', views.post_view, name='post_view'),
    path('post/search/', views.post_search, name='post_search'),

    path('new/post/', views.new_post, name='post_create'),
    path('new/comment/<str:uid>/', views.new_comment, name="create_comment"),


    path('b/list/', views.badge_list, name='badge_list'),
    path('t/list/', views.tags_list, name='tags_list'),

    path('b/view/<str:uid>/', views.badge_view, name='badge_view'),
    path('edit/post/<str:uid>/', views.edit_post, name='post_edit'),

    # Ajax calls
    path('ajax/users/search/', ajax.ajax_users_search, name='ajax_users_search'),
    path('ajax/digest/', ajax.ajax_digest, name='ajax_digest'),
    path('ajax/vote/', ajax.ajax_vote, name='vote'),
    path('ajax/test/', ajax.ajax_test, name='ajax_test'),
    path('ajax/tags/search/', ajax.ajax_tags_search, name='ajax_tags_search'),
    path('ajax/subscribe/', ajax.ajax_subs, name='ajax_sub'),
    path('drag/and/drop/', ajax.drag_and_drop, name='drag_and_drop'),
    path('similar/posts/<str:uid>/', ajax.similar_posts, name='similar_posts'),
    path('ajax/digest/<str:uid>/', ajax.ajax_digest, name='ajax_digest'),
    path('ajax/edit/<str:uid>/', ajax.ajax_edit, name='ajax_edit'),
    path('ajax/create/', ajax.ajax_create, name='ajax_create'),
    path('inplace/form/', ajax.inplace_form, name='inplace_form'),
    path('ajax/search/', ajax.ajax_search, name='ajax_search'),

    path('similar/posts/<str:uid>/', ajax.similar_posts, name='similar_posts'),

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




