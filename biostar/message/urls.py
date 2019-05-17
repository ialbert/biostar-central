from . import views
from django.conf.urls import url, include
import biostar.accounts.urls as account_patterns
urlpatterns = [

    # Message urls
    url(r'^$', views.index, name='index'),

    url(r'^inbox/$', views.message_list, name='inbox'),
    url(r'^outbox/$', views.outbox, name='outbox'),
    url(r'^inbox/(?P<uid>[-\w]+)/$', views.inbox_view, name='inbox_view'),
    url(r'^outbox/(?P<uid>[-\w]+)/$', views.outbox_view, name='outbox_view'),
    url(r'^reply/(?P<uid>[-\w]+)/$', views.reply, name="reply"),
    url(r'^compose/$', views.message_compose, name='compose'),

    # Include the accounts urls
    url(r'^accounts/', include(account_patterns)),

    ]