from django.conf.urls import url


from . import views


urlpatterns = [


    url(r'^post/list/$', views.post_list, name='post_list'),
    url(r'^post/create/$', views.post_create, name='post_create'),
    url(r'^post/view/(?P<uid>[-\w]+)$', views.post_view, name='post_view'),
    url(r'^post/edit/(?P<uid>[-\w]+)$', views.edit_post, name='post_edit'),
    url(r'^post/comment/(?P<uid>[-\w]+)$', views.post_comment, name='post_comment'),


]



