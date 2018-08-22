from . import views
from django.conf.urls import url

urlpatterns = [

    # Message urls
    url(r'^inbox/$', views.inbox_list, name='inbox'),
    url(r'^outbox/$', views.outbox_list, name='outbox'),
    url(r'^view/(?P<uid>[-\w]+)/$', views.message_view, name='message_view'),
    url(r'^compose/$', views.message_view, name='compose'),

    ]