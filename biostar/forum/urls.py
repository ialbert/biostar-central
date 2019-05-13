from . import views
from django.conf.urls import url

urlpatterns = [

    # Main entry. Post listing.
    url(r'^$', views.post_list, name='post_list'),

    url(r'^p/(?P<uid>[-\w]+)/$', views.post_view, name='post_view'),

    url(r'^b/list/$', views.badge_list, name='badge_list'),

    url(r'^b/view/(?P<uid>[-\w]+)/$', views.badge_view, name='badge_view'),

    url(r'^create/$', views.post_create, name='post_create'),
    url(r'^sub/(?P<uid>[-\w]+)/$', views.subs_action, name='subs_action'),
    url(r'^edit/(?P<uid>[-\w]+)/$', views.edit_post, name='post_edit'),
    url(r'^comment/$', views.comment, name='post_comment'),
    url(r"^comment/form/(?P<uid>[-\w]+)/$", views.comment_form, name="comment_form"),
    url(r'^vote/$', views.ajax_vote, name='vote'),
    url(r'^quick/feed/$', views.feed_post, name='feed_post'),

    #url(r'^tags/list/$', views.tags_list, name='tags_list'),
    url(r'^moderate/(?P<uid>[-\w]+)/$', views.post_moderate, name="post_moderate"),

    # Community urls
    url(r'^community/list/$', views.community_list, name='community_list'),

]






