from django.contrib import admin
from django.urls import path, include, re_path
import debug_toolbar
from django.conf import settings
from django.conf.urls.static import static
import biostar.accounts.urls as account_patterns
from biostar.recipes import views, api, ajax

urlpatterns = [
    # The django generated admin site.

    path('', views.index, name="index"),

    # Site
    path(r'site/admin/', views.site_admin, name='site_admin'),
    path(r'site/bin/', views.recycle_bin, name='recycle_bin'),

    # Ajax calls
    path(r'ajax/check/job/<str:uid>/', ajax.check_job, name='ajax_check_job'),
    #path(r'inplace/recipe/form/<str:uid>/', ajax.ajax_inplace_form, name='ajax_inplace_form'),
    path(r'add/recipe/fields/', ajax.add_to_interface, name="add_recipe_fields"),


    path(r'add/vars/', ajax.add_variables, name="add_vars"),
    path(r'preview/template/', ajax.preview_template, name="preview_template"),
    path(r'preview/json/', ajax.preview_json, name="preview_json"),

    path(r'toggle/delete/', ajax.toggle_delete, name="toggle_delete"),
    path(r'copy/object/', ajax.copy_object, name="copy_object"),
    path(r'manage/access/', ajax.manage_access, name="manage_access"),

    # Project
    path(r'project/users/<str:uid>/', views.project_users, name='project_users'),
    path(r'project/create/', views.project_create, name='project_create'),

    path(r'project/list/public/', views.project_list_public, name='project_list_public'),
    path(r'project/list/private/', views.project_list_private, name='project_list_private'),

    # This should not be needed/
    path(r'project/list/', views.project_list_public, name='project_list'),

    path(r'project/view/<str:uid>/', views.project_info, name='project_view'),

    path(r'<str:label>/view/', views.project_viewing, name='project_viewing'),

    path(r'project/edit/<str:uid>/', views.project_edit, name='project_edit'),
    path(r'project/info/<str:uid>/', views.project_info, name='project_info'),

    path(r'project/delete/<str:uid>/', views.project_delete, name='project_delete'),
    re_path(r'project/share/(?P<token>[-\w]+)/', views.project_share, name='project_share'),

    # Data
    path(r'<str:label>/data/', views.data_listing, name='data_listing'),
    path(r'data/list/<str:uid>/', views.data_list, name='data_list'),
    path(r'data/view/<str:uid>/', views.data_view, name='data_view'),
    path(r'data/edit/<str:uid>/', views.data_edit, name='data_edit'),
    path(r'data/upload/<str:uid>/', views.data_upload, name='data_upload'),
    re_path(r'^data/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.data_serve, name='data_serve'),
    path(r'data/paste/<str:uid>/', views.data_paste, name='data_paste'),
    path(r'data/delete/<str:uid>/', views.data_delete, name='data_delete'),

    # Recipes
    path(r'project/<str:label>/recipes/', views.recipe_listing, name='recipe_listing'),
    path(r'recipe/list/<str:uid>/', views.recipe_list, name='recipe_list'),

    path(r'latest/recipes/', views.latest_recipes, name='latest_recipes'),


    path(r'recipe/view/<str:uid>/', views.recipe_view, name='recipe_view'),
    path(r'recipe/run/<str:uid>/', views.recipe_run, name='recipe_run'),

    # Returns a rendered html fragment.
    path(r'get/part/<str:name>/<int:id>/', views.get_part, name='get_part'),

    path(r'ajax/recipe/edit/<int:id>/', ajax.ajax_edit, name='ajax_recipe_edit'),

    # Renders an HTML form field base on the TOML input.
    path(r'ajax/field/render/', ajax.field_render, name='ajax_field_render'),



    path(r'recipe/paste/<str:uid>/', views.recipe_paste, name='recipe_paste'),
    path(r'recipe/delete/<str:uid>/', views.recipe_delete, name='recipe_delete'),
    path(r'recipe/code/download/<str:uid>/', views.recipe_code_download, name='recipe_download'),
    path(r'recipe/create/<str:uid>/', views.recipe_create, name='recipe_create'),

    # File listings
    re_path(r'file/list/(?P<path>.+)$', views.import_files, name='file_list'),
    path(r'root/list/', views.import_files, name='root_list'),
    path(r'file/copy/', ajax.file_copy, name='file_copy'),

    # Actions
    path(r'action/clear/<str:uid>/', views.clear_clipboard, name='clear_clipboard'),
    path(r"search/", views.search_bar, name='search'),

    # Jobs
    path(r'<str:label>/results/', views.job_listing, name='job_listing'),
    path(r'job/list/<str:uid>/', views.job_list, name='job_list'),
    path(r'job/view/<str:uid>/', views.job_view, name='job_view'),
    path(r'job/edit/<str:uid>/', views.job_edit, name='job_edit'),
    re_path(r'^job/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.job_serve, name='job_serve'),
    path(r'job/delete/<str:uid>/', views.job_delete, name='job_delete'),
    path(r'job/rerun/<str:uid>/', views.job_rerun, name='job_rerun'),


    # Api calls
    path(r'api/list/', api.api_list, name='api_list'),

    path(r'api/recipe/json/<str:uid>/', api.recipe_json, name='recipe_api_json'),
    path(r'api/recipe/template/<str:uid>/', api.recipe_template, name='recipe_api_template'),
    path(r'api/recipe/image/<str:uid>/', api.recipe_image, name='recipe_api_image'),
    path(r'api/project/<str:uid>/', api.project_info, name='project_api_info'),
    path(r'api/project/image/<str:uid>/', api.project_image, name='project_api_image'),

    # Copy and paste actions
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

    urlpatterns += [
          path('__debug__/', include(debug_toolbar.urls)),
    ]


