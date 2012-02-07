from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns


# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('main.server',

    # index page
    (r'^$', 'views.index'),
    
    # show
    (r'^show/(?P<target>\w+)/$', 'views.index'),

    # static pages
    (r'^about/$','action.about'),
    (r'^tools/$', direct_to_template, {'template': 'tools.html'}),
    
    (r'^todo/$', direct_to_template, {'template': 'todo.html'}),
    (r'^feeds/$', direct_to_template, {'template': 'feeds.html'}),

    # moderation handlers
    (r'^cleanup/$', 'action.cleanup'),

    # tags
    (r'^tag/list/$', 'views.tag_list'),

    # badges
    (r'^badge/list/$', 'views.badge_list'),
    (r'^badge/show/(?P<bid>\d+)/$', 'action.badge_show'),
    
    # members
    (r'^user/list/$', 'views.user_list'),
    (r'^user/show/(?P<uid>\d+)/$', 'views.user_profile'),
    (r'^user/edit/(?P<uid>\d+)/$', 'action.user_edit'),

    # returns a preview page
    (r'^preview/$', 'views.preview'),
    
    # revisions    
    (r'^revisions/list/(?P<pid>\d+)/$', 'views.revision_list'),

    
    # questions
    (r'^question/unanswered/$', 'views.question_unanswered'),
    (r'^question/tagged/(?P<tag_name>[\w]+)/$', 'views.question_tagged'),
    
    # shows all the posts
    (r'^post/list/$', 'views.post_list'),
    (r'^post/list/(?P<word>[a-z\-]+)/$', 'views.post_list_filter'),

    (r'^post/list/(?P<uid>\d+)/$', 'views.post_list'),
    (r'^post/list/(?P<uid>\d+)/(?P<word>[a-z\-]+)/$', 'views.post_list_filter'),
  
    # post handlers with or withouth a slug
    (r'^post/show/(?P<pid>\d+)/$', 'views.post_show'),
    (r'^post/show/(?P<pid>\d+)/([-\w]+)/$', 'views.post_show'),
    
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/$','views.post_edit'),
    
    # handles new questions
    (r'^new/question/$','views.new_question'),
    
    # handles new answers
    (r'^new/answer/(?P<parentid>\d+)/$','views.new_answer'),
   
    # submits a new comment
    (r'^new/comment/(?P<parentid>\d+)/$','views.new_comment'),
    
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/(?P<parent>\d+)/$','views.post_edit'),
   
    # voting handler
    (r'^vote/$', 'views.vote'),
    
    # moderation handlers
    (r'^moderate/post/(?P<pid>\d+)/(?P<action>[a-z\-]+)/$', 'views.moderate_post'),
    (r'^moderate/user/(?P<uid>\d+)/(?P<action>[a-z\-]+)/$', 'views.moderate_user'),

    # destroys a post
    (r'^destroy/post/(?P<pid>\d+)/$', 'action.destroy_post'),
    

    # lists all moderator actions
    (r'^modlog/list/$', 'action.modlog_list'),
    
    # clear all notifications
    (r'^note/clear/(?P<uid>\d+)/$','action.note_clear'),
    
)

from server.feeds import LatestEntriesFeed, LatestNewsFeed

urlpatterns += patterns('',
    
    # RSS feeds
    (r'^feeds/latest/$', LatestEntriesFeed() ),
    (r'^feeds/messages/(?P<uuid>[a-z0-9]+)/$', LatestNewsFeed() ),

    # openid authentication
    (r'^openid/', include('django_openid_auth.urls')),
    (r'^logout/$', 'django.contrib.auth.views.logout',  {'next_page':'/'}),

    # Enable the admin:
    (r'^admin/', include(admin.site.urls)),
    (r'^admin/doc/', include('django.contrib.admindocs.urls')),
)

urlpatterns += staticfiles_urlpatterns()

 