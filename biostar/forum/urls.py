from . import views
from django.conf.urls import url

urlpatterns = [

    # Post urls
    url(r'^$', views.list_view, name='post_list'),
    url(r'^view/(?P<uid>[-\w]+)/$', views.post_view, name='post_view'),

    url(r'^create/$', views.post_create, name='post_create'),
    url(r'^sub/(?P<uid>[-\w]+)/$', views.subs_action, name='subs_action'),
    url(r'^edit/(?P<uid>[-\w]+)/$', views.edit_post, name='post_edit'),
    url(r'^comment/$', views.ajax_comment, name='post_comment'),
    url(r'^vote/(?P<uid>[-\w]+)/$', views.update_vote, name='update_vote'),
    url(r'^tags/list/$', views.tags_list, name='tags_list'),

    # Community urls
    url(r'^community/list/$', views.community_list, name='community_list'),

    # Message urls
    url(r'^messages/list/$', views.message_list, name='message_list'),
    url(r'^messages/view/(?P<uid>[-\w]+)/$', views.message_view, name='message_view'),

]






