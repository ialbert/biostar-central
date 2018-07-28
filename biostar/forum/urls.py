from . import views
from django.conf.urls import url
from django.conf import settings

msg_urls = [

    # Message urls
    url(r'^messages/list/$', views.message_list, name='message_list'),
    url(r'^messages/view/(?P<uid>[-\w]+)/$', views.message_view, name='message_view'),
]

urlpatterns = [

    # Post urls
    url(r'^$', views.list_view, name='post_list'),
    url(r'^view/(?P<uid>[-\w]+)/$', views.post_view, name='post_view'),

    url(r'^list/(?P<topic>[-\w]+)/$', views.list_by_topic, name='post_list_topic'),
    url(r'^create/$', views.post_create, name='post_create'),
    url(r'^sub/(?P<uid>[-\w]+)/$', views.subs_action, name='subs_action'),
    url(r'^edit/(?P<uid>[-\w]+)/$', views.edit_post, name='post_edit'),
    url(r'^comment/$', views.ajax_comment, name='post_comment'),
    url(r'^vote/(?P<uid>[-\w]+)/$', views.update_vote, name='update_vote'),

    # Community urls
    url(r'^community/list/$', views.community_list, name='community_list'),

]

# Add missing urls so templates do not break
if settings.ONLY_FORUM_URLS:
    urlpatterns += [

        url(r'^$', views.not_implemented, name='index'),
        url(r'^project/list/$', views.not_implemented, name='project_list'),
        url(r'^project/view/(?P<uid>[-\w]+)/$', views.not_implemented, name='project_view'),
        url(r'^data/view/(?P<uid>[-\w]+)/$', views.not_implemented, name='data_view'),
        url(r'^recipe/view/(?P<uid>[-\w]+)/$', views.not_implemented, name='recipe_view'),
        url(r'^job/view/(?P<uid>[-\w]+)/$', views.not_implemented, name='job_view'),
        url(r'^action/moderate/$', views.not_implemented, name='recipe_mod'),
        url(r'^site/bin/$', views.not_implemented, name='recycle_bin'),
        url(r'^site/admin/$', views.not_implemented, name='site_admin'),


    ]






