from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biostar3.forum.post_views import PostList, UserList, SearchResults, PostView

urlpatterns = patterns('',

    url(r'^$', PostList.as_view(), name='home'),

    # Renders search results.
    url(r'^search/$', SearchResults.as_view(), name='search'),

    # Post details.
    url(r'^p/(?P<pk>\d+)/$', PostView.as_view(), name="post_view"),

    # url(r'^blog/', include('blog.urls')),

    # The list of users.
    url(r'^user/list/$', UserList.as_view(), name="user_list"),

)
