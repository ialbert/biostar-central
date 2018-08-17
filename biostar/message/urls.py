from . import views
from django.conf.urls import url

urlpatterns = [

    # Message urls
    url(r'^list/$', views.message_list, name='message_list'),
    url(r'^view/(?P<uid>[-\w]+)/$', views.message_view, name='message_view'),
    url(r'^compose/$', views.message_view, name='compose'),

    ]