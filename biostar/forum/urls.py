from django.conf.urls import url


from . import views


urlpatterns = [

    url(r'^post/new/$', views.new_post, name='new_post'),
    url(r'^post/edit/(?P<uid>[-\w]+)$', views.edit_post, name='edit_post'),
    url(r'^post/answer/(?P<puid>[-\w]+)$', views.new_answer, name='new_answer')


]



