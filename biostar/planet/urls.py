
from django.urls import path
from biostar.planet import views

urlpatterns = [

    # Get the reset/ urls
    path('', views.blog_list, name="blog_list"),

]
