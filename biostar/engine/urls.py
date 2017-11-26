from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name="index"),

    url(r'^info/$', views.info, name="info"),

    # Engine specific admin site.
    url(r'^site/admin/', views.site_admin, name='site_admin'),

    url(r'^project/users/(?P<id>\d+)/$', views.project_users, name='project_users'),
    url(r'^project/create/$', views.project_create, name='project_create'),
    url(r'^project/list/$', views.project_list, name='project_list'),
    url(r'^project/view/(?P<id>\d+)/$', views.project_view, name='project_view'),
    url(r'^project/edit/(?P<id>\d+)/$', views.project_edit, name='project_edit'),

    url(r'^data/list/(?P<id>\d+)/$', views.data_list, name='data_list'),
    url(r'^data/view/(?P<id>\d+)/$', views.data_view, name='data_view'),
    url(r'^data/edit/(?P<id>\d+)/$', views.data_edit, name='data_edit'),
    url(r'^data/create/(?P<id>\d+)/$', views.data_upload, name='data_upload'),

    url(r'^analysis/list/(?P<id>\d+)/$', views.analysis_list, name='analysis_list'),
    url(r'^analysis/view/(?P<id>\d+)$', views.analysis_view, name='analysis_view'),
    url(r'^analysis/run/(?P<id>\d+)/$', views.analysis_run, name='analysis_run'),
    url(r'^analysis/recipe/(?P<id>\d+)$', views.analysis_recipe, name='analysis_recipe'),
    url(r'^recipe/edit/(?P<id>\d+)$', views.recipe_edit, name='recipe_edit'),
    url(r'^analysis/edit/(?P<id>\d+)$', views.analysis_edit, name='analysis_edit'),
    url(r'^analysis/copy/(?P<id>\d+)$', views.analysis_copy, name='analysis_copy'),

    url(r'^job/list/(?P<id>\d+)/$', views.job_list, name='job_list'),
    url(r'^job/view/(?P<id>\d+)/$', views.job_view, name='job_view'),
    url(r'^job/view/result/(?P<id>\d+)/$', views.job_result_view, name='job_result_view'),
    url(r'^job/view/files/(?P<id>\d+)/$', views.job_files_list, name='job_files_entry'),
    url(r'^job/view/files/(?P<id>\d+)/(?P<path>.+)/$', views.job_files_list, name='job_files_list'),

]

