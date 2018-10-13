from django.conf.urls import url

from . import views

urlpatterns = [

    # Get the reset/ urls
    url(r'^search/experiments$', views.search, name='exp_search')
    ]