from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns


# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',

    # index page
    (r'^$', 'main.server.views.index'),

    # static pages
    (r'^about/$','main.server.action.about'),
    (r'^tools/$', direct_to_template, {'template': 'tools.html'}),

    # tags
    (r'^tag/list/$', 'main.server.views.tag_list'),

    # badges
    (r'^badge/list/$', 'main.server.views.badge_list'),
    
    # search page
    (r'^search/$', 'main.server.views.search'),

    # members
    (r'^user/list/$', 'main.server.views.user_list'),
    (r'^user/show/(?P<uid>\d+)/$', 'main.server.views.user_profile'),
    (r'^user/edit/(?P<uid>\d+)/$', 'main.server.action.user_edit'),

    # returns a preview page
    (r'^preview/$', 'main.server.views.preview'),
    
    # revisions    
    (r'^revisions/(?P<pid>\d+)/list/$', 'main.server.views.revision_list'),

    # questions
    (r'^question/unanswered/$', 'main.server.views.question_unanswered'),
    (r'^question/list/$', 'main.server.views.question_list'),
    (r'^question/tagged/(?P<tag_name>[a-z\-]+)/$', 'main.server.views.question_tagged'),
    
    # post handlers with or withouth a slug
    (r'^post/show/(?P<pid>\d+)/$', 'main.server.views.post_show'),
    (r'^post/show/(?P<pid>\d+)/([-\w]+)/$', 'main.server.views.post_show'),
    
   
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/$','main.server.views.post_edit'),
    
    # handles new questions
    (r'^new/question/$','main.server.views.new_question'),
    
    # handles new answers
    (r'^new/answer/(?P<parentid>\d+)/$','main.server.views.new_answer'),
   
    # submits a new comment
    (r'^new/comment/(?P<parentid>\d+)/$','main.server.views.new_comment'),
    
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/(?P<parent>\d+)/$','main.server.views.post_edit'),
   
    # voting handler
    (r'^vote/$', 'main.server.views.vote'),
    (r'^moderate/$', 'main.server.views.moderate'),
    (r'^modlog/list/$', 'main.server.action.modlog_list'),
    (r'^note/clear/(?P<uid>\d+)/$','main.server.action.note_clear'),
    

    # openid authentication
    (r'^openid/', include('django_openid_auth.urls')),
    (r'^logout/$', 'django.contrib.auth.views.logout',  {'next_page':'/'}),

    # Enable the admin:
    (r'^admin/', include(admin.site.urls)),
    (r'^admin/doc/', include('django.contrib.admindocs.urls')),
)

urlpatterns += staticfiles_urlpatterns()

 