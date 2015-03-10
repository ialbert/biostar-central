from __future__ import absolute_import, division, print_function, unicode_literals
from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings

urlpatterns = patterns('',

    # The main forum.
    url(r'^', include('biostar3.forum.urls')),

    # The authentication backend.
    (r'^accounts/', include('allauth.urls')),

    # Admin url pattern.
    url(r'^admin/', include(admin.site.urls)),

)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
