from django.contrib import admin
from django.urls import path, include, re_path
import debug_toolbar
from django.conf import settings
from django.conf.urls.static import static
from biostar.accounts.urls import account_patterns
from biostar.accounts.views import image_upload_view
from biostar.recipes import views, api, ajax

recipes_patterns = [

    path('', views.index, name="index"),

    # Site
    path(r'site/admin/', views.site_admin, name='site_admin'),
    path(r'site/bin/', views.recycle_bin, name='recycle_bin'),

    # Ajax calls
    path(r'ajax/check/job/<str:uid>/', ajax.check_job, name='ajax_check_job'),
    path(r'preview/json/', ajax.preview_json, name="preview_json"),
    path(r'toggle/delete/', ajax.toggle_delete, name="toggle_delete"),
    path(r'manage/access/', ajax.manage_access, name="manage_access"),
    path(r'recipe/drop/', ajax.drop_recipe, name="drop_recipe"),
    path(r'project/drop/', ajax.drop_project, name="drop_project"),

    # Ajax clipboard actions.
    path(r'clear/', ajax.ajax_clear_clipboard, name='clear_clipboard'),
    path(r'file/copy/', ajax.copy_file, name='copy_file'),
    path(r'copy/object/', ajax.copy_object, name="copy_object"),
    path(r'clipboard/', ajax.ajax_clipboard, name="ajax_clipboard"),
    path(r'paste/', ajax.ajax_paste, name='ajax_paste'),

    # Project
    path(r'project/users/<str:uid>/', views.project_users, name='project_users'),
    path(r'project/create/', views.project_create, name='project_create'),

    # Redirects users to public or private projects.
    path(r'project/list/', views.project_list, name='project_list'),


    path(r'project/view/<str:uid>/', views.project_info, name='project_view'),
    path(r'project/edit/<str:uid>/', views.project_edit, name='project_edit'),
    path(r'project/info/<str:uid>/', views.project_info, name='project_info'),
    path(r'project/delete/<str:uid>/', views.project_delete, name='project_delete'),
    re_path(r'project/share/(?P<token>[-\w]+)/', views.project_share, name='project_share'),

    # Data
    path(r'data/list/<str:uid>/', views.data_list, name='data_list'),
    path(r'data/view/<str:uid>/', views.data_view, name='data_view'),
    path(r'data/edit/<str:uid>/', views.data_edit, name='data_edit'),
    path(r'data/upload/<str:uid>/', views.data_upload, name='data_upload'),
    path(r'data/download/<str:uid>/', views.data_download, name='data_download'),
    path(r'data/delete/<str:uid>/', views.data_delete, name='data_delete'),
    re_path(r'^data/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.data_serve, name='data_serve'),

    # Recipes
    path(r'recipe/list/<str:uid>/', views.recipe_list, name='recipe_list'),
    path(r'latest/recipes/', views.latest_recipes, name='latest_recipes'),
    path(r'recipe/view/<str:uid>/', views.recipe_view, name='recipe_view'),
    path(r'recipe/run/<str:uid>/', views.recipe_run, name='recipe_run'),

    # Returns a rendered html fragment.
    path(r'get/part/<str:name>/<int:id>/', views.get_part, name='get_part'),
    path(r'ajax/recipe/edit/<int:id>/', ajax.ajax_edit, name='ajax_recipe_edit'),
    # Renders an HTML form field base on the TOML input.
    path(r'ajax/field/render/', ajax.field_render, name='ajax_field_render'),
    path(r'ajax/move/', ajax.ajax_move, name='ajax_move'),

    path(r'recipe/delete/<str:uid>/', views.recipe_delete, name='recipe_delete'),
    path(r'recipe/code/download/<str:uid>/<str:fname>', views.recipe_code_download, name='recipe_download'),
    path(r'recipe/create/<str:uid>/', views.recipe_create, name='recipe_create'),

    # File listings
    re_path(r'^file/list/(?P<path>.+)$', views.import_files, name='file_list'),
    path(r'root/list/', views.import_files, name='root_list'),
    # Actions
    path(r"search/", views.search_bar, name='search'),

    # Jobs
    path(r'job/list/<str:uid>/', views.job_list, name='job_list'),
    path(r'job/view/<str:uid>/', views.job_view, name='job_view'),
    path(r'job/edit/<str:uid>/', views.job_edit, name='job_edit'),
    re_path(r'^job/serve/(?P<uid>[-\w]+)/(?P<path>.+)$', views.job_serve, name='job_serve'),
    path(r'job/delete/<str:uid>/', views.job_delete, name='job_delete'),
    path(r'job/rerun/<str:uid>/', views.job_rerun, name='job_rerun'),

    # Api calls
    path(r'api/list/', api.api_list, name='api_list'),
    path(r'api/project/', api.project_api, name='project_api'),
    path(r'api/recipe/', api.recipe_api, name='recipe_api'),
    path(r'api/data/', api.data_api, name='data_api'),

    # Plugins
    path(r'render/plugin/', ajax.render_plugins, name='render_plugins'),

]


urlpatterns = [

    path(r'', include(recipes_patterns)),

    # Include the accounts urls
    path(r'accounts/', include(account_patterns)),

]

if settings.PAGEDOWN_IMAGE_UPLOAD_ENABLED:

    urlpatterns += [
        # Pagedown image upload url.
        path('pagedown/image-upload/', image_upload_view, name="pagedown-image-upload")
    ]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)

    urlpatterns += [
          path('__debug__/', include(debug_toolbar.urls)),
    ]

