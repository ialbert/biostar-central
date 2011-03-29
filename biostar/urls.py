from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',

    # index page
    (r'^$', 'biostar.server.views.index'),

    ('^about/$', direct_to_template, {'template': 'about.html'}),
    ('^newquestion/$','biostar.server.views.newquestion'),


    (r'^members/$', 'biostar.server.views.users'),
    (r'^member/(?P<uid>\d+)/$', 'biostar.server.views.user'),
    (r'^question/(?P<pid>\d+)/$', 'biostar.server.views.question'),
    (r'^newpost/$', 'biostar.server.views.newpost'),
    
    (r'^vote/$', 'biostar.server.views.vote'),

    # openid authentication
    (r'^openid/', include('django_openid_auth.urls')),
    (r'^logout/$', 'django.contrib.auth.views.logout',  {'next_page':'/'}),

    # Example:
    # (r'^biostar/', include('biostar.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # (r'^admin/', include(admin.site.urls)),

    (r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_DIR }),
)
