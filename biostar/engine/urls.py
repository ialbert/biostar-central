from django.conf.urls import url
from django.urls import path

from biostar.forum import views as forum_views
from . import views, api

urlpatterns = [
    url(r'^$', views.index, name="index"),

    # Site
    url(r'^site/admin/$', views.site_admin, name='site_admin'),
    url(r'^site/bin/$', views.recycle_bin, name='recycle_bin'),

    # Project
    url(r'^project/users/(?P<uid>[-\w]+)/$', views.project_users, name='project_users'),
    url(r'^project/create/$', views.project_create, name='project_create'),
    url(r'^project/list/$', views.project_list, name='project_list'),
    url(r'^project/view/(?P<uid>[-\w]+)/$', views.project_info, name='project_view'),
    url(r'^project/edit/(?P<uid>[-\w]+)/$', views.project_edit, name='project_edit'),
    url(r'^project/info/(?P<uid>[-\w]+)/$', views.project_info, name='project_info'),
    url(r'^project/list/private/$', views.project_list_private, name='project_list_private'),
    url(r'^project/list/public/$', views.project_list, name='project_list_public'),
    url(r'^project/delete/(?P<uid>[-\w]+)/$', views.project_delete, name='project_delete'),

    # Data
    url(r'^data/list/(?P<uid>[-\w]+)/$', views.data_list, name='data_list'),
    url(r'^data/view/(?P<uid>[-\w]+)/$', views.data_view, name='data_view'),
    url(r'^data/edit/(?P<uid>[-\w]+)/$', views.data_edit, name='data_edit'),
    url(r'^data/upload/(?P<uid>[-\w]+)/$', views.data_upload, name='data_upload'),
    url(r'^data/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.data_serve, name='data_serve'),
    url(r'^data/paste/(?P<uid>[-\w]+)/$', views.data_paste, name='data_paste'),
    url(r'^data/delete/(?P<uid>[-\w]+)/$', views.data_delete, name='data_delete'),

    # Recipes
    url(r'^recipe/list/(?P<uid>[-\w]+)/$', views.recipe_list, name='recipe_list'),
    url(r'^recipe/view/(?P<uid>[-\w]+)/$', views.recipe_view, name='recipe_view'),
    url(r'^recipe/run/(?P<uid>[-\w]+)/$', views.recipe_run, name='recipe_run'),
    url(r'^recipe/edit/(?P<uid>[-\w]+)/$', views.recipe_edit, name='recipe_edit'),
    url(r'^recipe/code/view/(?P<uid>[-\w]+)/$', views.recipe_code_view, name='recipe_code_view'),
    url(r'^recipe/code/edit/(?P<uid>[-\w]+)/$', views.recipe_code_edit, name='recipe_code_edit'),
    url(r'^recipe/paste/(?P<uid>[-\w]+)/$', views.recipe_paste, name='recipe_paste'),
    url(r'^recipe/delete/(?P<uid>[-\w]+)/$', views.recipe_delete, name='recipe_delete'),
    url(r'^recipe/code/download/(?P<uid>[-\w]+)/$', views.recipe_code_download, name='recipe_download'),

    # Actions
    url(r'^action/clear/(?P<uid>[-\w]+)/$', views.clear_clipboard, name='clear_clipboard'),
    url(r'^action/subscribe/(?P<uid>[-\w]+)/$', views.discussion_subs, name='discussion_subs'),

    url(r"^search/$", views.search_bar, name='search'),


    # Jobs
    url(r'^job/list/(?P<uid>[-\w]+)/$', views.job_list, name='job_list'),
    url(r'^job/view/(?P<uid>[-\w]+)/$', views.job_view, name='job_view'),
    url(r'^job/edit/(?P<uid>[-\w]+)/$', views.job_edit, name='job_edit'),
    url(r'^job/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.job_serve, name='job_serve'),
    url(r'^job/delete/(?P<uid>[-\w]+)/$', views.job_delete, name='job_delete'),

    # Api calls
    url(r'^recipe/api/list/$', api.recipe_api_list, name='recipe_api_list'),
    url(r'^project/api/list/$', api.project_api_list, name='project_api_list'),
    url(r'^api/recipe/(?P<uid>[-\w]+)/json/$', api.recipe_json, name='recipe_api_json'),
    url(r'^api/recipe/(?P<uid>[-\w]+)/template/$', api.recipe_template, name='recipe_api_template'),
    url(r'^api/recipe/(?P<uid>[-\w]+)/image/$', api.recipe_image, name='recipe_api_image'),

    # Discussions
    url(r'^discussion/list/(?P<uid>[-\w]+)/$', views.discussion_list, name='discussion_list'),
    url(r'^discussion/create/(?P<uid>[-\w]+)/$', views.discussion_create, name='discussion_create'),
    url(r'^discussion/view/(?P<uid>[-\w]+)/$', views.discussion_view, name='discussion_view'),
    url(r'^comment/$', forum_views.comment, name='discussion_comment'),

    # Ajax calls
    #url(r'^data/copy/$', views.ajax_data_copy, name='data_copy'),
    #url(r'^result/copy/$', views.ajax_job_copy, name='job_copy'),
    #url(r'^recipe/copy/$', views.ajax_recipe_copy, name='recipe_copy'),

    url(r'^data/copy/(?P<uid>[-\w]+)/$', views.data_copy, name='data_copy'),
    url(r'^result/copy/(?P<uid>[-\w]+)/$', views.job_copy, name='job_copy'),
    url(r'^recipe/copy/(?P<uid>[-\w]+)/$', views.recipe_copy, name='recipe_copy'),
    url(r'^data/file/copy/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.data_file_copy, name='data_file_copy'),
    url(r'^job/file/copy/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.job_file_copy, name='job_file_copy'),
    url(r'^file/paste/(?P<uid>[-\w]+)/$', views.file_paste, name='file_paste'),

]

