from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',

    # index page
    (r'^$', 'biostar.server.views.index'),

    # static pages
    ('^about/$', direct_to_template, {'template': 'about.html'}),

    # tags
    (r'^tag/list/$', 'biostar.server.views.tag_list'),

    # badges
    (r'^badge/list/$', 'biostar.server.views.badge_list'),
    
    # search page
    (r'^search/$', 'biostar.server.views.search'),

    # members
    (r'^member/list/$', 'biostar.server.views.user_list'),
    (r'^member/show/(?P<uid>\d+)/$', 'biostar.server.views.user_profile'),

    # returns a preview page
    (r'^preview/$', 'biostar.server.views.preview'),
    
    # revisions
    
    (r'^revisions/(?P<pid>\d+)/list/$', 'biostar.server.views.revision_list'),

    
    # questions
    (r'^question/unanswered/$', 'biostar.server.views.question_unanswered'),
    (r'^question/list/$', 'biostar.server.views.question_list'),
    (r'^question/tagged/(?P<tag_name>[a-z\-]+)/$', 'biostar.server.views.question_tagged'),
    (r'^question/show/(?P<pid>\d+)/$', 'biostar.server.views.question_show'),
    (r'^question/edit/(?P<pid>\d+)/$', 'biostar.server.views.question_edit'),
    (r'^question/new/$','biostar.server.views.question_edit'),
    
    # answers
    (r'^answer/edit/(?P<qid>\d+)/(?P<aid>\d+)/$', 'biostar.server.views.answer_edit'),
    (r'^answer/new/(?P<qid>\d+)/$', 'biostar.server.views.answer_edit'),
    
    # comment handlers
    (r'^comment/new/(?P<pid>\d+)/$', 'biostar.server.views.comment_add'),

    # voting handler
    (r'^vote/$', 'biostar.server.views.vote'),

    # openid authentication
    (r'^openid/', include('django_openid_auth.urls')),
    (r'^logout/$', 'django.contrib.auth.views.logout',  {'next_page':'/'}),

    # Enable the admin:
    (r'^admin/', include(admin.site.urls)),
    (r'^admin/doc/', include('django.contrib.admindocs.urls')),

   
)

if settings.ADMIN_PASSWORD_OVERRIDE:
    urlpatterns += patterns('', 
        # admin login override
        (r'^admin/password/override/$', 'biostar.server.views.admin_password_override'),
    )

if settings.DEBUG:
    urlpatterns += patterns('', 
        # static content
        (r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_DIR }),
    )
 