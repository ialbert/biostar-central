from django.conf.urls import url, patterns, include
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from main.server import const

from server.views_refactored import MessageView

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

handler500 = 'main.server.action.url500'

urlpatterns = patterns('main.server',

    # main index page
    url(r'^$', 'views.index', name='index'),
    
    url(r'^search/$', 'search.main', name="search"),
    url(r'^more/like/(?P<pid>\d+)/$', 'search.more', name="more"),

    # shows new messages
    url(r'^show/messages/$', MessageView.as_view(), name=MessageView.url),


    # show by content type
    url(r'^show/(?P<tab>\w+)/$', 'views.index', name="show"),
    
    # show tagged posts
    url(r'^show/tag/(?P<tag_name>.+)/$', 'views.show_tag', name="show-tag"),
    
    # show posts by user
    url(r'^show/user/(?P<uid>\d+)/$', 'views.show_user', name="show-user"),
    url(r'^show/user/(?P<uid>\d+)/(?P<post_type>\w+)/$', 'views.show_user', name="show-user-content"),
    url(r'^linkout/(?P<pid>\d+)/$', 'views.linkout', name="linkout"),
    
    # urls for the navigation bar
    url(r'^tag/list/$', 'views.tag_list', name="tag-list"),
    url(r'^user/list/$', 'views.user_list', name="user-list"),
    url(r'^badge/list/$', 'views.badge_list', name="badge-list"),
    
    
    # user edit page
    url(r'^user/edit/(?P<uid>\d+)/$', 'action.user_edit', name="user-edit"),
    
    # new style user profile
    url(r'^u/(?P<uid>\d+)/$', 'views.user_profile', name="user-profile"),
    url(r'^u/(?P<uid>\d+)/(?P<tab>\w+)/$', 'views.user_profile', name="user-profile-tab"),

    # old style user profile
    url(r'^user/profile/(?P<uid>\d+)/$', 'views.user_profile_redirect', name="user-profile-redirect"),
    url(r'^user/profile/(?P<uid>\d+)/(?P<tab>\w+)/$', 'views.user_profile_redirect', name="user-profile-tab-redirect"),
        
    # moderation handlers
    #(r'^cleanup/$', 'action.cleanup'),

    # badges
    url(r'^badge/show/(?P<bid>\d+)/$', 'action.badge_show', name="badge-show"),
    
    # revisions    
    url(r'^revision/show/(?P<pid>\d+)/$', 'views.revision_show', name="revision-show"),

    # new-style (short) post handlers
    url(r'^p/(?P<pid>\d+)/$', 'views.post_show', name="post-show"),

    # old style post show
    url(r'^post/show/(?P<pid>\d+)/$', 'views.post_show_redirect', name="post-show-redirect"),
    url(r'^post/show/(?P<pid>\d+)/([-\w]+)/$', 'views.post_show_redirect', name="post-show-slug-redirect"),
    url(r'^post/redirect/(?P<pid>\d+)/$', 'views.post_redirect', name="post-redirect"),
    
    # turned off reparenting for now
    #url(r'^post/reparent/(?P<pid>\d+)/$', 'action.post_reparent', name="post-reparent"),

    # editing an existing post/answer/comment
    url(r'^post/edit/(?P<pid>\d+)/$','views.post_edit', name="post-edit"),

    # private messaging
    url(r'^pm/(?P<uid>\d+)/$','action.private_message', name="private-message"),

    # moderation views
    url(r'^user/moderate/(?P<uid>\d+)/(?P<status>\w+)/$','action.user_moderate', name="user-moderate"),    
    url(r'^post/moderate/(?P<pid>\d+)/(?P<status>\w+)/$','action.post_moderate', name="post-moderate"),
    url(r'^merge/$','action.request_merge', name="request-merge"),
    url(r'^approve_merge/(?P<master_id>\d+)/(?P<remove_id>\d+)/$','action.approve_merge', name="approve-merge"),
    url(r'^request/info/(?P<pid>\d+)/$','pages.request_info', name="request-info"),
    
    # handles new post
    url(r'^new/post/$','views.new_post', name="new-post"),
    url(r'^new/answer/(?P<pid>\d+)/$','views.new_answer', name="new-answer"),
    url(r'^new/comment/(?P<pid>\d+)/$','views.new_comment', name="new-comment"),
    
    # static pages
    url(r'^about/$','pages.about', name='about'),
    url(r'^rss/$','pages.rss', name='rss'),
    url(r'^faq/$','pages.faq', name='faq'),
    url(r'^beta/$','pages.beta', name='beta'),
    url(r'^google/$','pages.google', name='google'),
    url(r'^testpage/$','pages.testpage', name='testpage'),

    # help pages go here
    url(r'^help/x/$','pages.help_external'),

    # lists all moderator actions
    url(r'^modlog/list/$', 'views.modlog_list', name="modlog-list"),
  
    # destroys a post
    url(r'^comment/delete/(?P<pid>\d+)/$', 'ajax.comment_delete', name="comment-delete"),
      
    # voting handler
    url(r'^vote/$', 'ajax.vote', name="vote"),
    url(r'^tagcomplete/$', 'ajax.tagcomplete', name="tagcomplete"),


    # clear all notifications
    url(r'^note/clear/(?P<uid>\d+)/$','action.note_clear', name="note-clear"),
   
    # redirecting to new post
    url(r'^questions/(?P<pid>\d+)/$','action.redirect_post', name="redirect-short"),
    url(r'^questions/(?P<pid>\d+)/([-\w]+)/$','action.redirect_post', name="redirect-post"),
    url(r'^questions/tagged/(?P<tag>.+)/$','action.redirect_tag', name="redirect-tag"),




    # the main handler for the external authentication
    url(r'^x/$','action.external_handler', name="external-handler"),

    #url(r'^x/post/$','action.external_post', name="external-post"),

    # test login, used during debugging
    url(r'^test/login/(?P<uid>\d+)/(?P<token>[\w\d]+)/$','action.test_login', name="test-login"),
   
    # json api for stat generation
    url(r'^api/traffic/$', 'api.traffic', name='stats-traffic'),

    url(r'^api/user/(?P<uid>\d+)/$', 'api.user_info', name='api-user'),
    url(r'^api/post/(?P<pid>\d+)/$', 'api.post_info', name='api-post'),
    url(r'^api/stats/(?P<days>\d+)/$', 'api.stats', name='api-stats'),

)



#
# Generic views
#
from django.views.generic.list import ListView
from django.views.generic import TemplateView
from main.server import models

urlpatterns += patterns('',
    url(r'^blog/list/$', ListView.as_view(
        queryset = models.Blog.objects.all().select_related('author__profile'),
        template_name='generic/blog.list.html'), name='blog-list'),
    
     # matching the robots.txt
    url(r'^robots\.txt$', TemplateView.as_view(template_name="robots.txt"), name='robots'),

)

#
# RSS Feeds 
#
from server.feeds import LatestEntriesFeed, NotificationFeed, MyTagsFeed, PostTypeFeed
from server.feeds import TagsFeed, PostFeed, UserFeed

urlpatterns += patterns('',
    
    # RSS feeds
    url(r'^feeds/latest/$', LatestEntriesFeed(), name='latest-feed' ),
    url(r'^feeds/messages/(?P<uuid>[a-z0-9]+)/$', NotificationFeed(), name='notification-feed' ),
    url(r'^feeds/mytags/(?P<uuid>[a-z0-9]+)/$', MyTagsFeed(), name='mytags-feed' ),
    url(r'^feeds/tag/(?P<text>[\w\-_\+]+)/$', TagsFeed(), name='tags-feed' ),
    url(r'^feeds/post/(?P<text>[\w\-_\+]+)/$', PostFeed(), name='post-feed' ),
    url(r'^feeds/user/(?P<text>[\w\-_\+]+)/$', UserFeed(), name='user-feed' ),
    url(r'^feeds/type/(?P<text>[\w\-_\+]+)/$', PostTypeFeed(), name='post-type-feed' ),
    
    # openid authentication
    url(r'^openid/', include('django_openid_auth.urls'), name='openid-login'),
    url(r'^logout/$', 'django.contrib.auth.views.logout',  {'next_page':'/'}, name='logout'),

    # Enable the admin:
    url(r'^admin/', include(admin.site.urls), name='admin'),
    url(r'^admin/doc/', include('django.contrib.admindocs.urls'), name='admin-docs'),
)

urlpatterns += staticfiles_urlpatterns()

 