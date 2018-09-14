from django.conf.urls import url
from django.urls import path

from biostar.forum import views as forum_views
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
    url(r'^data/paste/(?P<uid>[-\w]+)/$', views.data_paste, name='data_paste'),

    # Recipes
    url(r'^recipe/list/(?P<uid>[-\w]+)/$', views.recipe_list, name='recipe_list'),
    url(r'^recipe/view/(?P<uid>[-\w]+)/$', views.recipe_view, name='recipe_view'),
    url(r'^recipe/create/(?P<uid>[-\w]+)/$', views.recipe_create, name='recipe_create'),
    url(r'^recipe/run/(?P<uid>[-\w]+)/$', views.recipe_run, name='recipe_run'),
    url(r'^recipe/edit/(?P<uid>[-\w]+)/$', views.recipe_edit, name='recipe_edit'),
    url(r'^recipe/code/view/(?P<uid>[-\w]+)/$', views.recipe_code_view, name='recipe_code_view'),
    url(r'^recipe/code/edit/(?P<uid>[-\w]+)/$', views.recipe_code_edit, name='recipe_code_edit'),
    url(r'^recipe/diff/(?P<uid>[-\w]+)/$', views.recipe_diff, name='recipe_diff'),
    url(r'^recipe/paste/(?P<uid>[-\w]+)/$', views.recipe_paste, name='recipe_paste'),

    # Actions
    url(r'^action/clear/(?P<uid>[-\w]+)/$', views.clear_clipboard, name='clear_clipboard'),
    url(r'^action/toggle/(?P<uid>[-\w]+)/(?P<obj_type>[-\w]+)/$', views.object_state_toggle, name='toggle_state'),
    url(r'^action/moderate/$', views.recipe_mod, name='recipe_mod'),
    url(r'^action/subscribe/(?P<uid>[-\w]+)/$', views.discussion_subs, name='discussion_subs'),

    # Jobs
    url(r'^job/list/(?P<uid>[-\w]+)/$', views.job_list, name='job_list'),
    url(r'^job/view/(?P<uid>[-\w]+)/$', views.job_view, name='job_view'),
    url(r'^job/edit/(?P<uid>[-\w]+)/$', views.job_edit, name='job_edit'),
    url(r'^job/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.job_serve, name='job_serve'),

    # Discussions
    url(r'^discussion/list/(?P<uid>[-\w]+)/$', views.discussion_list, name='discussion_list'),
    url(r'^discussion/create/(?P<uid>[-\w]+)/$', views.discussion_create, name='discussion_create'),
    url(r'^discussion/view/(?P<uid>[-\w]+)/$', views.discussion_view, name='discussion_view'),
    url(r'^discussion/comment/$', forum_views.comment, name='discussion_comment'),

    # Ajax calls
    url(r'^data/copy/$', views.ajax_data_copy, name='data_copy'),
    url(r'^result/copy/$', views.ajax_job_copy, name='job_copy'),
    url(r'^recipe/copy/$', views.ajax_recipe_copy, name='recipe_copy'),

    # Apis
    url(r'^api/recipe/(?P<uid>[-\w]+)/$', api.recipe_details, name='recipe_details'),


]

