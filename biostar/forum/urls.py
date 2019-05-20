
from django.conf import settings
from django.conf.urls.static import static
from django.conf.urls import url, include

from biostar.forum import views
import biostar.message.urls as msg_patterns
import biostar.accounts.urls as account_patterns

urlpatterns = [

    # Main entry. Post listing.
    url(r'^$', views.latest, name='post_list'),
    url(r'^myvotes/$', views.myvotes, name='myvotes'),
    url(r'^bookmarks/$', views.bookmarks, name='bookmarks'),
    url(r'^following/$', views.following, name='following'),
    url(r'^myposts/$', views.myposts, name='myposts'),

    url(r'^p/(?P<uid>[-\w]+)/$', views.post_view, name='post_view'),

    url(r'^b/list/$', views.badge_list, name='badge_list'),

    url(r'^b/view/(?P<uid>[-\w]+)/$', views.badge_view, name='badge_view'),

    url(r'^create/$', views.post_create, name='post_create'),
    url(r'^sub/(?P<uid>[-\w]+)/$', views.subs_action, name='subs_action'),
    url(r'^edit/(?P<uid>[-\w]+)/$', views.edit_post, name='post_edit'),
    url(r'^comment/$', views.comment, name='post_comment'),
    url(r"^comment/form/(?P<uid>[-\w]+)/$", views.comment_form, name="comment_form"),
    url(r'^vote/$', views.ajax_vote, name='vote'),

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





