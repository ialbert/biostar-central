from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name="index"),

    url(r'^docs/(?P<name>[-\w]+)/$', views.docs, name='docs'),

    # Engine specific admin site.
    url(r'^site/admin/', views.site_admin, name='site_admin'),

    url(r'^project/users/(?P<uid>[-\w]+)/$', views.project_users, name='project_users'),
    url(r'^project/create/$', views.project_create, name='project_create'),
    url(r'^project/list/$', views.project_list, name='project_list'),
    url(r'^project/view/(?P<uid>[-\w]+)/$', views.project_view, name='project_view'),
    url(r'^project/edit/(?P<uid>[-\w]+)/$', views.project_edit, name='project_edit'),
    url(r'^project/datatypes/(?P<uid>[-\w]+)/$', views.project_types, name='project_types'),

    url(r'^data/list/(?P<uid>[-\w]+)/$', views.data_list, name='data_list'),
    url(r'^data/view/(?P<id>\d+)/$', views.data_view, name='data_view'),
    url(r'^data/edit/(?P<id>\d+)/$', views.data_edit, name='data_edit'),
    url(r'^data/upload/(?P<uid>[-\w]+)/$', views.data_upload, name='data_upload'),

    # Recipe URLS
    url(r'^recipe/list/(?P<uid>[-\w]+)/$', views.recipe_list, name='recipe_list'),
    url(r'^recipe/view/(?P<id>\d+)/$', views.recipe_view, name='recipe_view'),
    url(r'^recipe/create/(?P<uid>[-\w]+)/$', views.recipe_create, name='recipe_create'),
    url(r'^recipe/run/(?P<id>\d+)/$', views.recipe_run, name='analysis_run'),
    url(r'^recipe/edit/(?P<id>\d+)/$', views.recipe_edit, name='recipe_edit'),
    url(r'^recipe/code/(?P<id>\d+)/$', views.recipe_code, name='recipe_code'),

    url(r'^job/list/(?P<uid>[-\w]+)/$', views.job_list, name='job_list'),
    url(r'^job/view/(?P<id>\d+)/$', views.job_view, name='job_view'),
    url(r'^job/edit/(?P<id>\d+)/$', views.job_edit, name='job_edit'),
    url(r'^job/view/result/(?P<id>\d+)/$', views.job_result_view, name='job_result_view'),
    url(r'^job/view/files/(?P<id>\d+)/$', views.job_files_list, name='job_files_entry'),
    url(r'^job/view/files/(?P<id>\d+)/(?P<path>.+)/$', views.job_files_list, name='job_files_list'),

    # Block the media dir from being exposed
    url(r'^media/$', views.block_media_url, name="block"),
    url(r'^media/(?P<path1>[-\w]+)/$', views.block_media_url, name="block"),
    url(r'^media/(?P<path1>[-\w]+)/(?P<path2>[-\w]+)/$', views.block_media_url, name="block")

]

