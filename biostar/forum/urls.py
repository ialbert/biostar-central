from django.conf.urls import url


from . import views


urlpatterns = [


    url(r'^post/list/$', views.list_view, name='post_list'),
    url(r'^messages/list/$', views.message_list, name='message_list'),
    url(r'^post/create/$', views.post_create, name='post_create'),
    url(r'^post/subs/action/(?P<uid>[-\w]+)$', views.subs_action, name='subs_action'),
    url(r'^post/view/(?P<uid>[-\w]+)$', views.post_view, name='post_view'),
    url(r'^post/edit/(?P<uid>[-\w]+)$', views.edit_post, name='post_edit'),
    url(r'^post/comment/(?P<uid>[-\w]+)$', views.post_comment, name='post_comment'),
    url(r'^post/vote/(?P<uid>[-\w]+)$', views.update_vote, name='update_vote'),


]



