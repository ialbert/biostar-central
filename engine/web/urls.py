from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name="home"),
    url(r'^blog/$', views.blog, name="blog"),
    url(r'^data/$', views.data, name="data"),
    url(r'^forum/$', views.forum, name="forum"),
]
