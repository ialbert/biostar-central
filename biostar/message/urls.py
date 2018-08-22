from . import views
from django.conf.urls import url

urlpatterns = [

    # Message urls
    url(r'^inbox/$', views.inbox_list, name='inbox'),
    url(r'^outbox/$', views.outbox_list, name='outbox'),
    url(r'^inbox/(?P<uid>[-\w]+)/$', views.inbox_view, name='inbox_view'),
    url(r'^outbox/(?P<uid>[-\w]+)/$', views.outbox_view, name='outbox_view'),
    url(r'^compose/$', views.message_view, name='compose'),

    ]