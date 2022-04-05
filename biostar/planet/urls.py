
from django.urls import path, include
from biostar.planet import views
from django.contrib import admin


planet_patterns = [
    path('', views.blog_list, name="blog_list"),
    path(r'view/<int:id>/', views.blog_view, name="blog_view"),
    path(r'bump/<int:id>/', views.blog_bump, name="blog_bump"),
]

urlpatterns = [

    # Get the reset/ urls
    path(r'', include(planet_patterns)),
    # Add admin urls.
    path('admin/', admin.site.urls),

]
