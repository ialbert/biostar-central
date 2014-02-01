from django.conf.urls import patterns, include, url
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

from biostar.server import views

urlpatterns = patterns('',
    # Homepage.
    url(r'^$', views.Index.as_view(), name="home"),

    # The list of users.
    url(r'^user/list/$', views.UserList.as_view(), name="userlist"),

    # User profile.
    url(r'^u/(?P<pk>\d+)/$', views.UserProfile.as_view(), name="profile"),

    # Social login pages.
    (r'^accounts/', include('allauth.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
)


# Adding the flatpages.
urlpatterns += patterns('django.contrib.flatpages.views',
    url(r'^info/about/$', 'flatpage', {'url': '/about/'}, name='about'),
    url(r'^info/faq/$', 'flatpage', {'url': '/faq/'}, name='faq'),
    url(r'^info/help/$', 'flatpage', {'url': '/about/'}, name='help'),
    url(r'^info/policy/$', 'flatpage', {'url': '/policy/'}, name='policy'),
)

# This is used only for the debug toolbar
if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
        url(r'^__debug__/', include(debug_toolbar.urls)),
    )