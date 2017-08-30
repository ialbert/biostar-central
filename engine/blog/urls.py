from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index),
    url(r'^blog/$', views.blog),
    url(r'^data/$', views.data),
    url(r'^forum/$', views.forum),
]
