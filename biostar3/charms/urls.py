from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.conf import settings
from . import views

urlpatterns = patterns('',
    # Biostar charms.
    url(r'^$', views.charms, name="charm_view"),
)
