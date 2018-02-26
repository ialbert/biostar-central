from django.conf.urls import url
from django.urls import path

from . import views

urlpatterns = [
    url(r'^$', views.index, name="index"),

    # url(r'^docs/(?P<name>[-\w]+)/$', views.docs, name='docs'),

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
    url(r'^data/copy/(?P<uid>[-\w]+)/$', views.data_copy, name='data_copy'),
    url(r'^data/paste/(?P<uid>[-\w]+)/$', views.data_paste, name='data_paste'),

    # Data file browsing.
    url(r'^data/navigate/(?P<uid>[-\w]+)/$', views.data_navigate, name='data_navigate'),
    url(r'^data/browse/(?P<uid>[-\w]+)/$', views.data_browse, name='data_entry'),
    url(r'^data/browse/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.data_browse, name='data_browser'),

    # Data file serve.
    url(r'^data/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.data_serve, name='data_serve'),

    # Recipes
    url(r'^recipe/list/(?P<uid>[-\w]+)/$', views.recipe_list, name='recipe_list'),
    url(r'^recipe/view/(?P<uid>[-\w]+)/$', views.recipe_view, name='recipe_view'),
    url(r'^recipe/create/(?P<uid>[-\w]+)/$', views.recipe_create, name='recipe_create'),
    url(r'^recipe/run/(?P<uid>[-\w]+)/$', views.recipe_run, name='recipe_run'),
    url(r'^recipe/edit/(?P<uid>[-\w]+)/$', views.recipe_edit, name='recipe_edit'),
    url(r'^recipe/code/(?P<uid>[-\w]+)/$', views.recipe_code, name='recipe_code'),
    url(r'^recipe/copy/(?P<uid>[-\w]+)/$', views.recipe_copy, name='recipe_copy'),
    url(r'^recipe/paste/(?P<uid>[-\w]+)/$', views.recipe_paste, name='recipe_paste'),
    url(r'^recipe/paste/(?P<uid>[-\w]+)/$', views.recipe_paste, name='recipe_paste'),
    url(r'^recipe/diff/(?P<uid>[-\w]+)/$', views.recipe_diff, name='recipe_diff'),

    # Actions
    url(r'^action/clear/(?P<uid>[-\w]+)/(?P<url>.+)/(?P<board>.+)/$', views.clear_clipboard, name='clear_clipboard'),
    url(r'^action/toggle/(?P<uid>[-\w]+)/(?P<obj_type>[-\w]+)/$', views.object_state_toggle, name='toggle_state'),
    url(r'^action/paste/(?P<uid>[-\w]+)/$', views.files_paste, name='files_paste'),

    # Jobs
    url(r'^job/list/(?P<uid>[-\w]+)/$', views.job_list, name='job_list'),
    url(r'^job/view/(?P<uid>[-\w]+)/$', views.job_view, name='job_view'),
    url(r'^job/edit/(?P<uid>[-\w]+)/$', views.job_edit, name='job_edit'),

    # Job file browsing.
    url(r'^job/browse/(?P<uid>[-\w]+)/$', views.job_browse, name='job_entry'),
    url(r'^job/browse/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.job_browse, name='job_browser'),

    # Job file serve.
    url(r'^job/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.job_serve, name='job_serve')

]

