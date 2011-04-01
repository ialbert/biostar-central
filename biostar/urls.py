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
    
    (r'^members/$', 'biostar.server.views.user_list'),
    (r'^members/show/(?P<uid>\d+)/$', 'biostar.server.views.user_profile'),
    
    # question handlers
    (r'^question/list/$', 'biostar.server.views.question_list'),
    (r'^question/show/(?P<pid>\d+)/$', 'biostar.server.views.question_show'),
    (r'^question/edit/(?P<pid>\d+)/$', 'biostar.server.views.question_edit'),
    (r'^question/new/$','biostar.server.views.question_edit'),
    
    # answer handlers, questionid, andswerid
    (r'^answer/edit/(?P<qid>\d+)/(?P<aid>\d+)/$', 'biostar.server.views.answer_edit'),
    (r'^answer/new/(?P<qid>\d+)/$', 'biostar.server.views.answer_edit'),
    
    # comment handlers
    (r'^comment/new/(?P<pid>\d+)/$', 'biostar.server.views.comment_add'),


    # voting handler
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
