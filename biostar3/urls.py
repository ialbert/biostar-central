from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum.views import PostList, UserList

urlpatterns = patterns('',

    url(r'^$', PostList.as_view(), name='home'),

    # url(r'^blog/', include('blog.urls')),

    # The list of users.

    url(r'^user/list/$', UserList.as_view(), name="user-list"),

    url(r'^admin/', include(admin.site.urls)),
)
