from django.conf.urls import url, include
from django.conf import settings
from . import views

urlpatterns = [

    # Get the reset/ urls
    #url('^', include('django.contrib.auth.urls')),

    url(r'^$', views.index, name='emailer-index'),

    ]
