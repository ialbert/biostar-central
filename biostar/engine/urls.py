from django.contrib import admin
from django.urls import path, include, re_path

from django.conf import settings
from django.conf.urls.static import static
import biostar.accounts.urls as account_patterns
from biostar.engine import views, api

urlpatterns = [
    # The django generated admin site.

    path('', views.index, name="index"),

    # Site
    path(r'site/admin/', views.site_admin, name='site_admin'),
    path(r'site/bin/', views.recycle_bin, name='recycle_bin'),

    # Project
    path(r'project/users/<str:uid>/', views.project_users, name='project_users'),
    path(r'project/create/', views.project_create, name='project_create'),
    path(r'project/list/', views.project_list, name='project_list'),
    path(r'project/view/<str:uid>/', views.project_info, name='project_view'),
    path(r'project/edit/<str:uid>/', views.project_edit, name='project_edit'),
    path(r'project/info/<str:uid>/', views.project_info, name='project_info'),
    path(r'project/list/private/', views.project_list_private, name='project_list_private'),
    path(r'project/list/public/', views.project_list_public, name='project_list_public'),
    path(r'project/list/', views.project_list, name='project_list'),
    path(r'project/delete/<str:uid>/', views.project_delete, name='project_delete'),

    # Data
    path(r'data/list/<str:uid>/', views.data_list, name='data_list'),
    path(r'data/view/<str:uid>/', views.data_view, name='data_view'),
    path(r'data/edit/<str:uid>/', views.data_edit, name='data_edit'),
    path(r'data/upload/<str:uid>/', views.data_upload, name='data_upload'),
    re_path(r'^data/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.data_serve, name='data_serve'),
    path(r'data/paste/<str:uid>/', views.data_paste, name='data_paste'),
    path(r'data/delete/<str:uid>/', views.data_delete, name='data_delete'),

    # Recipes
    path(r'recipe/list/<str:uid>/', views.recipe_list, name='recipe_list'),
    path(r'recipe/view/<str:uid>/', views.recipe_view, name='recipe_view'),
    path(r'recipe/run/<str:uid>/', views.recipe_run, name='recipe_run'),
    path(r'recipe/edit/<str:uid>/', views.recipe_edit, name='recipe_edit'),
    #path(r'^recipe/code/view/(?P<uid>[-\w]+)/$', views.recipe_code_view, name='recipe_code_view'),
    path(r'recipe/code/edit/<str:uid>/', views.recipe_code_edit, name='recipe_code_edit'),
    path(r'recipe/paste/<str:uid>/', views.recipe_paste, name='recipe_paste'),
    path(r'recipe/delete/<str:uid>/', views.recipe_delete, name='recipe_delete'),
    path(r'recipe/code/download/<str:uid>/', views.recipe_code_download, name='recipe_download'),

    # Actions
    path(r'action/clear/<str:uid>/', views.clear_clipboard, name='clear_clipboard'),
    path(r"search/", views.search_bar, name='search'),

    # Jobs
    path(r'job/list/<str:uid>/', views.job_list, name='job_list'),
    path(r'job/view/<str:uid>/', views.job_view, name='job_view'),
    path(r'job/edit/<str:uid>/', views.job_edit, name='job_edit'),
    re_path(r'^job/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.job_serve, name='job_serve'),
    path(r'job/delete/<str:uid>/', views.job_delete, name='job_delete'),

    # Api calls
    path(r'api/list/', api.api_list, name='api_list'),

    path(r'api/recipe/json/<str:uid>/', api.recipe_json, name='recipe_api_json'),
    path(r'api/recipe/template/<str:uid>/', api.recipe_template, name='recipe_api_template'),
    path(r'api/recipe/image/<str:uid>/', api.recipe_image, name='recipe_api_image'),
    path(r'api/project/<str:uid>/', api.project_info, name='project_api_info'),
    path(r'api/project/image/<str:uid>/', api.project_image, name='project_api_image'),

    path(r'data/copy/<str:uid>/', views.data_copy, name='data_copy'),
    path(r'result/copy/<str:uid>/', views.job_copy, name='job_copy'),
    path(r'recipe/copy/<str:uid>/', views.recipe_copy, name='recipe_copy'),
    re_path(r'^data/file/copy/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.data_file_copy, name='data_file_copy'),
    re_path(r'^job/file/copy/(?P<uid>[-\w]+)/(?P<path>.+)/$', views.job_file_copy, name='job_file_copy'),
    path(r'file/paste/<str:uid>/', views.file_paste, name='file_paste'),

    # Include the accounts urls
    path(r'accounts/', include(account_patterns)),

]


if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)




