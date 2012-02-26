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
    (r'^show/tag/(?P<tag_name>[\w\-_]+)/$', 'views.show_tag'),
    
    # show posts by user
    (r'^show/user/(?P<uid>\d+)/$', 'views.show_user'),
    (r'^show/user/(?P<uid>\d+)/(?P<post_type>\w+)/$', 'views.show_user'),
    
    # urls for the navigation bar
    (r'^tag/list/$', 'views.tag_list'),
    (r'^user/list/$', 'views.user_list'),
    (r'^badge/list/$', 'views.badge_list'),
    
    
    # user edit page
    (r'^user/edit/(?P<uid>\d+)/$', 'action.user_edit'),
    
    # show user profile
    (r'^user/profile/(?P<uid>\d+)/$', 'views.user_profile'),
    (r'^user/profile/(?P<uid>\d+)/(?P<tab>\w+)/$', 'views.user_profile'),
    (r'^user/moderate/(?P<uid>\d+)/(?P<status>\w+)/$','action.user_moderate'),
    
    #
    # old handlers
    #
    (r'^tools/$', direct_to_template, {'template': 'tools.html'}),
    (r'^todo/$', direct_to_template, {'template': 'todo.html'}),
   
    # moderation handlers
    (r'^cleanup/$', 'action.cleanup'),

   
    # badges
    (r'^badge/show/(?P<bid>\d+)/$', 'action.badge_show'),
    
    # revisions    
    (r'^revision/show/(?P<pid>\d+)/$', 'views.revision_show'),

    
    # post handlers with or withouth a slug
    (r'^post/show/(?P<pid>\d+)/$', 'views.post_show'),
    (r'^post/show/(?P<pid>\d+)/([-\w]+)/$', 'views.post_show'),
    
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/$','views.post_edit'),
    (r'^post/moderate/(?P<pid>\d+)/(?P<status>\w+)/$','action.post_moderate'),
    
    # handles new post
    (r'^new/post/$','views.new_post'),
    (r'^new/answer/(?P<pid>\d+)/$','views.new_answer'),
    (r'^new/comment/(?P<pid>\d+)/$','views.new_comment'),
    

    # static pages
    (r'^about/$','pages.about'),
    (r'^rss/$','pages.rss'),
    (r'^faq/$','pages.faq'),
    (r'^beta/$','pages.beta'),
    
    # lists all moderator actions
    (r'^modlog/list/$', 'views.modlog_list'),
  

    # ----------------

    # destroys a post
    (r'^post/destroy/(?P<pid>\d+)/$', 'ajax.post_destroy'),
    
    # ajax handlers
    
    # returns a preview page
    (r'^preview/$', 'ajax.preview'),
    
    # voting handler
    (r'^vote/$', 'ajax.vote'),
    
     (r'^moderate/user/(?P<uid>\d+)/(?P<action>[a-z\-]+)/$', 'ajax.moderate_user'),

    

   
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

 