from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name="index"),
    url(r'^docs/(?P<name>[-\w]+)/$', views.docs, name='docs'),

    # Engine specific admin site.
    url(r'^site/admin/$', views.site_admin, name='site_admin'),

    url(r'^recycle/bin/$', views.recycle_bin, name='recycle_bin'),

    url(r'^project/users/(?P<uid>[-\w]+)/$', views.project_users, name='project_users'),
    url(r'^project/create/$', views.project_create, name='project_create'),
    url(r'^project/list/$', views.project_list, name='project_list'),
    url(r'^project/view/(?P<uid>[-\w]+)/$', views.project_view, name='project_view'),
    url(r'^project/edit/(?P<uid>[-\w]+)/$', views.project_edit, name='project_edit'),

    url(r'^data/list/(?P<uid>[-\w]+)/$', views.data_list, name='data_list'),
    url(r'^data/view/(?P<uid>[-\w]+)/$', views.data_view, name='data_view'),
    url(r'^data/edit/(?P<uid>[-\w]+)/$', views.data_edit, name='data_edit'),
    url(r'^data/upload/(?P<uid>[-\w]+)/$', views.data_upload, name='data_upload'),
    url(r'^data/view/files/(?P<uid>[-\w]+)/$', views.data_files_list, name='data_files_entry'),
    url(r'^data/view/files/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.data_files_list, name='data_files_list'),
    url(r'^data/copy/(?P<uid>[-\w]+)/$', views.data_copy, name='data_copy'),
    url(r'^data/paste/(?P<uid>[-\w]+)/$', views.data_paste, name='data_paste'),
    url(r'^files/paste/(?P<uid>[-\w]+)/$', views.files_paste, name='files_paste'),
    url(r'^data/navigate/(?P<uid>[-\w]+)/$', views.data_nav, name='data_nav'),
    url(r'^data/file/serve/(?P<uid>[-\w]+)/(?P<file_path>.+)/$', views.data_file_serve, name='data_file_serve'),
    url(r'^data/delete/(?P<uid>[-\w]+)/(?P<state>[-\w]+)/$', views.data_state_change, name='data_delete'),
    url(r'^data/restore/(?P<uid>[-\w]+)/(?P<state>[-\w]+)/$', views.data_state_change, name='data_restore'),

    # Recipe URLS
    url(r'^recipe/list/(?P<uid>[-\w]+)/$', views.recipe_list, name='recipe_list'),
    url(r'^recipe/view/(?P<uid>[-\w]+)/$', views.recipe_view, name='recipe_view'),
    url(r'^recipe/create/(?P<uid>[-\w]+)/$', views.recipe_create, name='recipe_create'),
    url(r'^recipe/run/(?P<uid>[-\w]+)/$', views.recipe_run, name='recipe_run'),
    url(r'^recipe/edit/(?P<uid>[-\w]+)/$', views.recipe_edit, name='recipe_edit'),
    url(r'^recipe/code/(?P<uid>[-\w]+)/$', views.recipe_code, name='recipe_code'),
    url(r'^recipe/copy/(?P<uid>[-\w]+)/$', views.recipe_copy, name='recipe_copy'),
    url(r'^recipe/paste/(?P<uid>[-\w]+)/$', views.recipe_paste, name='recipe_paste'),

    url(r'^clear/clipboard/(?P<uid>[-\w]+)/(?P<url>.+)/(?P<board>.+)/$', views.clear_clipboard, name='clear_clipboard'),

    url(r'^job/list/(?P<uid>[-\w]+)/$', views.job_list, name='job_list'),
    url(r'^job/view/(?P<uid>[-\w]+)/$', views.job_view, name='job_view'),
    url(r'^job/edit/(?P<uid>[-\w]+)/$', views.job_edit, name='job_edit'),
    url(r'^job/delete/(?P<uid>[-\w]+)/(?P<state>[-\w]+)/$', views.job_state_change, name='job_delete'),
    url(r'^job/restore/(?P<uid>[-\w]+)/(?P<state>[-\w]+)/$', views.job_state_change, name='job_restore'),
    url(r'^job/view/files/(?P<uid>[-\w]+)/$', views.job_files_list, name='job_files_entry'),
    url(r'^job/view/files/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.job_files_list, name='job_files_list'),
    url(r'^job/file/serve/(?P<uid>[-\w]+)/(?P<file_path>.+)/$', views.job_file_serve, name='job_file_serve')

]

