

from django.conf.urls import url, include
from . import views

urlpatterns = [

    # Get the reset/ urls
    #url('^', include('django.contrib.auth.urls')),

    url(r'^password/reset/$', views.password_reset, name='password_reset'),
    url(r'^password/reset/done/$', views.password_reset_done, name='password_reset_done'),

    url(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        views.pass_reset_confirm, name='password_reset_confirm'),
    url(r'^reset/done/$', views.password_reset_complete, name='password_reset_complete'),

    url(r'^login/$', views.user_login, name="login"),
    url(r'^login/$', views.user_login, name="login"),
    url(r'^signup/$', views.user_signup, name="signup"),
    url(r'^profile/(?P<id>\d+)/$', views.profile, name="profile"),
    url(r'^edit/profile/(?P<id>\d+)/$', views.edit_profile, name='edit_profile'),
    url(r'^logout/$', views.user_logout, name="logout")

]

