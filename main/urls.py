from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns


# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('main.server',

    # main index page
    (r'^$', 'views.index'),
    
    # show by content type
    (r'^show/(?P<tab>\w+)/$', 'views.index'),
    
    # show tagged posts
    (r'^show/tag/(?P<tag_name>[\w]+)/$', 'views.show_tag'),
    
    # show posts by user
    (r'^show/user/(?P<uid>\d+)/$', 'views.show_user'),
    (r'^show/user/(?P<uid>\d+)/(?P<post_type>\w+)/$', 'views.show_user'),
    
    
    # urls for the navigation bar
    (r'^tag/list/$', 'views.tag_list'),
    (r'^user/list/$', 'views.user_list'),
    (r'^badge/list/$', 'views.badge_list'),
    (r'^about/$','action.about'),
    (r'^rss/$','action.rss'),
    (r'^faq/$','action.faq'),
    
    
    # user edit page
    (r'^user/edit/(?P<uid>\d+)/$', 'action.user_edit'),
    
    # show user profile
    (r'^user/profile/(?P<uid>\d+)/$', 'views.user_profile'),
    (r'^user/profile/(?P<uid>\d+)/(?P<tab>\w+)/$', 'views.user_profile'),
    #
    # old handlers
    #
    (r'^tools/$', direct_to_template, {'template': 'tools.html'}),
    
    (r'^todo/$', direct_to_template, {'template': 'todo.html'}),
   
    # moderation handlers
    (r'^cleanup/$', 'action.cleanup'),

    

    # badges
   
    (r'^badge/show/(?P<bid>\d+)/$', 'action.badge_show'),
    
    # members
    
   
    
    # revisions    
    (r'^revisions/list/(?P<pid>\d+)/$', 'views.revision_list'),

    
    # post handlers with or withouth a slug
    (r'^post/show/(?P<pid>\d+)/$', 'views.post_show'),
    (r'^post/show/(?P<pid>\d+)/([-\w]+)/$', 'views.post_show'),
    
    
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/$','views.post_edit'),
    
    # handles new post
    (r'^new/post/$','views.new_question'),
    
    
    # handles new questions
    (r'^new/question/$','views.new_question'),
    
    # handles new answers
    (r'^new/answer/(?P<parentid>\d+)/$','views.new_answer'),
   
    # submits a new comment
    (r'^new/comment/(?P<parentid>\d+)/$','views.new_comment'),
    
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/(?P<parent>\d+)/$','views.post_edit'),
   
    # ajax handlers
    
    # returns a preview page
    (r'^preview/$', 'ajax.preview'),
    
    # voting handler
    (r'^vote/$', 'ajax.vote'),
    
    # moderation handlers
    (r'^moderate/post/(?P<pid>\d+)/(?P<action>[a-z\-]+)/$', 'ajax.moderate_post'),
    (r'^moderate/user/(?P<uid>\d+)/(?P<action>[a-z\-]+)/$', 'ajax.moderate_user'),

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

 