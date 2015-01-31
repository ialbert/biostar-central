from __future__ import (unicode_literals, division, absolute_import)
from django.conf.urls import patterns, include, url

from . import views

urlpatterns = patterns('',
    url(r'^$', views.TorrentList.as_view(), name='index'),
    url(r'^peers/$', views.PeerList.as_view(), name='peer-list'),
    url(r'^download/(?P<pk>\d+)/$', views.download_torrent, name='download'),
    url(r'^view/peer/(?P<pk>\d+)/$', views.PeerDetail.as_view(), name='view-peer'),
    url(r'^torrent/view/(?P<pk>\d+)/$', views.TorrentDetail.as_view(), name='view-torrent'),
    url(r'^torrent/delete/(?P<pk>\d+)/$', views.TorrentDelete.as_view(), name='torrent-delete'),
    url(r'^announce$', views.announce, name='announce'),
    url(r'^scrape$', views.scrape, name='scrape'),
)