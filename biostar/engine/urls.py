from django.conf.urls import url
from django.urls import path

from . import views
from . import api

urlpatterns = [
    url(r'^$', views.index, name="index"),

    # Site
    url(r'^site/admin/$', views.site_admin, name='site_admin'),
    url(r'^site/bin/$', views.recycle_bin, name='recycle_bin'),

    # Project
    url(r'^project/users/(?P<uid>[-\w]+)/$', views.project_users, name='project_users'),
    url(r'^project/create/$', views.project_create, name='project_create'),
    url(r'^project/list/$', views.project_list, name='project_list'),
    url(r'^project/view/(?P<uid>[-\w]+)/$', views.project_view, name='project_view'),
    url(r'^project/edit/(?P<uid>[-\w]+)/$', views.project_edit, name='project_edit'),

    # Data
    url(r'^data/list/(?P<uid>[-\w]+)/$', views.data_list, name='data_list'),
    url(r'^data/view/(?P<uid>[-\w]+)/$', views.data_view, name='data_view'),
    url(r'^data/edit/(?P<uid>[-\w]+)/$', views.data_edit, name='data_edit'),
    url(r'^data/upload/(?P<uid>[-\w]+)/$', views.data_upload, name='data_upload'),
    url(r'^data/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.data_serve, name='data_serve'),

    # Recipes
    url(r'^recipe/list/(?P<uid>[-\w]+)/$', views.recipe_list, name='recipe_list'),
    url(r'^recipe/view/(?P<uid>[-\w]+)/$', views.recipe_view, name='recipe_view'),
    url(r'^recipe/create/(?P<uid>[-\w]+)/$', views.recipe_create, name='recipe_create'),
    url(r'^recipe/run/(?P<uid>[-\w]+)/$', views.recipe_run, name='recipe_run'),
    url(r'^recipe/edit/(?P<uid>[-\w]+)/$', views.recipe_edit, name='recipe_edit'),
    url(r'^recipe/code/(?P<uid>[-\w]+)/$', views.recipe_code, name='recipe_code'),
    url(r'^recipe/diff/(?P<uid>[-\w]+)/$', views.recipe_diff, name='recipe_diff'),

    # Actions
    url(r'^action/clear/(?P<uid>[-\w]+)/(?P<url>.+)/(?P<board>.+)/$', views.clear_clipboard, name='clear_clipboard'),
    url(r'^action/toggle/(?P<uid>[-\w]+)/(?P<obj_type>[-\w]+)/$', views.object_state_toggle, name='toggle_state'),
    url(r'^action/moderate/$', views.recipe_mod, name='recipe_mod'),

    # Jobs
    url(r'^job/list/(?P<uid>[-\w]+)/$', views.job_list, name='job_list'),
    url(r'^job/view/(?P<uid>[-\w]+)/$', views.job_view, name='job_view'),
    url(r'^job/edit/(?P<uid>[-\w]+)/$', views.job_edit, name='job_edit'),
    url(r'^job/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.job_serve, name='job_serve'),

    # Discussions
    url(r'^discussion/list/(?P<uid>[-\w]+)/$', views.discussion_list, name='discussion_list'),
    url(r'^discussion/create/(?P<uid>[-\w]+)/$', views.discussion_create, name='discussion_create'),

    # Apis
    url(r'^api/recipe/(?P<uid>[-\w]+)/$', api.recipe_details, name='recipe_details'),


]

