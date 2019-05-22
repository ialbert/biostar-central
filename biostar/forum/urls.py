
from django.conf import settings
from django.conf.urls.static import static
from django.conf.urls import url, include

from biostar.forum import views
from biostar.forum import ajax
import biostar.message.urls as msg_patterns
import biostar.accounts.urls as account_patterns

urlpatterns = [

    # Main entry. Post listing.
    url(r'^$', views.latest, name='post_list'),
    url(r'^votes/$', views.myvotes, name='myvotes'),
    url(r'^bookmarks/$', views.bookmarks, name='bookmarks'),
    url(r'^following/$', views.following, name='following'),
    url(r'^myposts/$', views.myposts, name='myposts'),

    url(r'^p/(?P<uid>[-\w]+)/$', views.post_view, name='post_view'),


    url(r'^create/post/$', views.post_create, name='post_create'),
    url(r'^create/answer/(?P<uid>[-\w]+)/$', views.post_answer, name='post_answer'),
    url(r"^create/comment/(?P<uid>[-\w]+)/$", views.create_comment, name="create_comment"),


    url(r'^b/list/$', views.badge_list, name='badge_list'),

    url(r'^b/view/(?P<uid>[-\w]+)/$', views.badge_view, name='badge_view'),


    url(r'^sub/(?P<uid>[-\w]+)/$', views.subs_action, name='subs_action'),
    url(r'^edit/(?P<uid>[-\w]+)/$', views.edit_post, name='post_edit'),

    url(r'^ajax/vote/$', ajax.ajax_vote, name='vote'),

    url(r'^ajax/test/$', ajax.ajax_test, name='ajax_test'),

    #url(r'^tags/list/$', views.tags_list, name='tags_list'),
    url(r'^moderate/(?P<uid>[-\w]+)/$', views.post_moderate, name="post_moderate"),

    # Community urls
    url(r'^community/$', views.community_list, name='community_list'),

    # Include messages urls
    url(r'^message/', include(msg_patterns)),

    # Include the accounts urls
    url(r'^accounts/', include(account_patterns)),

]


if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)





