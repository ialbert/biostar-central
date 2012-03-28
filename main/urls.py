from django.conf.urls.defaults import *
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('main.server',

    # main index page
    (r'^$', 'views.index'),
    
    (r'^search/$', 'search.main'),
    (r'^more/like/(?P<pid>\d+)/$', 'search.more'),
    
    # show by content type
    (r'^show/(?P<tab>\w+)/$', 'views.index'),
    
    # show tagged posts
    (r'^show/tag/(?P<tag_name>[\w\-_\+]+)/$', 'views.show_tag'),
    
    # show posts by user
    (r'^show/user/(?P<uid>\d+)/$', 'views.show_user'),
    (r'^show/user/(?P<uid>\d+)/(?P<post_type>\w+)/$', 'views.show_user'),
    (r'^show/blog/(?P<pid>\d+)/$', 'views.blog_redirect'),
    
    # urls for the navigation bar
    (r'^tag/list/$', 'views.tag_list'),
    (r'^user/list/$', 'views.user_list'),
    (r'^badge/list/$', 'views.badge_list'),
    
    
    # user edit page
    (r'^user/edit/(?P<uid>\d+)/$', 'action.user_edit'),
    
    # show user profile
    (r'^user/profile/(?P<uid>\d+)/$', 'views.user_profile'),
    (r'^user/profile/(?P<uid>\d+)/(?P<tab>\w+)/$', 'views.user_profile'),
    
    
    # moderation handlers
    #(r'^cleanup/$', 'action.cleanup'),

   
    # badges
    (r'^badge/show/(?P<bid>\d+)/$', 'action.badge_show'),
    
    # revisions    
    (r'^revision/show/(?P<pid>\d+)/$', 'views.revision_show'),

    # post handlers with or withouth a slug
    (r'^post/show/(?P<pid>\d+)/$', 'views.post_show'),
    (r'^post/show/(?P<pid>\d+)/([-\w]+)/$', 'views.post_show'),
    (r'^post/redirect/(?P<pid>\d+)/$', 'views.post_redirect'),
    
    # editing an existing post/answer/comment
    (r'^post/edit/(?P<pid>\d+)/$','views.post_edit'),
    
    # moderation views
    (r'^user/moderate/(?P<uid>\d+)/(?P<status>\w+)/$','action.user_moderate'),    
    (r'^post/moderate/(?P<pid>\d+)/(?P<status>\w+)/$','action.post_moderate'),
    (r'^merge/$','action.request_merge'),
    (r'^approve_merge/(?P<master_id>\d+)/(?P<remove_id>\d+)/$','action.approve_merge'),
    
    # handles new post
    (r'^new/post/$','views.new_post'),
    (r'^new/answer/(?P<pid>\d+)/$','views.new_answer'),
    (r'^new/comment/(?P<pid>\d+)/$','views.new_comment'),
    
    
    
    # static pages
    (r'^about/$','pages.about'),
    (r'^rss/$','pages.rss'),
    (r'^faq/$','pages.faq'),
    (r'^beta/$','pages.beta'),
    (r'^google/$','pages.google'),
    
    # lists all moderator actions
    (r'^modlog/list/$', 'views.modlog_list'),
  

    # ----------------
    
    # destroys a post
    (r'^comment/delete/(?P<pid>\d+)/$', 'ajax.comment_delete'),
    
    # ajax handlers
    
    # voting handler
    (r'^vote/$', 'ajax.vote'),
    
       
    # clear all notifications
    (r'^note/clear/(?P<uid>\d+)/$','action.note_clear'),
   
)

#
# Generic views
#
from django.views.generic.list import ListView
from main.server import models

urlpatterns += patterns('',
    (r'^blog/list/$', ListView.as_view(
        queryset = models.Blog.objects.all().select_related('author__profile'),
        template_name='generic/blog.list.html')),    
)

#
# RSS Feeds 
#
from server.feeds import LatestEntriesFeed, NotificationFeed, MyTagsFeed
from server.feeds import TagsFeed, PostFeed, UserFeed

urlpatterns += patterns('',
    
    # RSS feeds
    (r'^feeds/latest/$', LatestEntriesFeed() ),
    (r'^feeds/messages/(?P<uuid>[a-z0-9]+)/$', NotificationFeed() ),
    (r'^feeds/mytags/(?P<uuid>[a-z0-9]+)/$', MyTagsFeed() ),
    (r'^feeds/tag/(?P<text>[\w\-_\+]+)/$', TagsFeed() ),
    (r'^feeds/post/(?P<text>[\w\-_\+]+)/$', PostFeed() ),
    (r'^feeds/user/(?P<text>[\w\-_\+]+)/$', UserFeed() ),

    # openid authentication
    (r'^openid/', include('django_openid_auth.urls')),
    (r'^logout/$', 'django.contrib.auth.views.logout',  {'next_page':'/'}),

    # Enable the admin:
    (r'^admin/', include(admin.site.urls)),
    (r'^admin/doc/', include('django.contrib.admindocs.urls')),
)

urlpatterns += staticfiles_urlpatterns()

 