from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

from biostar.server import views

urlpatterns = patterns('',
    # Homepage.
    url(r'^$', views.IndexView.as_view(), name="home"),

    # The list of users.
    url(r'^user/list/$', views.UserView.as_view(), name="userlist"),

    # Social login pages.
    (r'^accounts/', include('allauth.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),

)
