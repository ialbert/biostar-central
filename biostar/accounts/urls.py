
from django.conf.urls import url, include
from django.conf import settings
from . import views

urlpatterns = [

    # Get the reset/ urls
    #url('^', include('django.contrib.auth.urls')),

    url(r'^password/reset/$', views.password_reset, name='password_reset'),
    url(r'^password/reset/done/$', views.password_reset_done, name='password_reset_done'),

    url(r'^verify/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        views.email_verify_account, name='email_verify_account'),

    url(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        views.pass_reset_confirm, name='password_reset_confirm'),

    url(r'^verify/$', views.send_email_verify, name="send_email_verify"),

    url(r'^reset/done/$', views.password_reset_complete, name='password_reset_complete'),
    url(r'^moderate/(?P<uid>[-\w]+)/$', views.user_moderate, name="user_moderate"),
    url(r'^login/$', views.user_login, name="login"),
    url(r'^signup/$', views.user_signup, name="signup"),
    url(r'^profile/(?P<uid>[-\w]+)/$', views.user_profile, name="user_profile"),

    url(r'^edit/profile/$', views.edit_profile, name='edit_profile'),
    url(r'^toggle/notify/$', views.toggle_notify, name='toggle_notify'),
    url(r'^logout/$', views.user_logout, name="logout"),

    url(r'^debug/user$', views.debug_user, name="debug_user"),

    # External url login
    url(r'^external/$', views.external_login, name="external"),

    # Used for 3rd party logins.
    url("^social/", include('allauth.urls')),

]

if settings.ALLOW_SELF_MODERATE:
    # Allow users to toggle their moderation state
    urlpatterns += [url(r'^toggle/moderate/$', views.toggle_moderate, name="toggle_moderate")]


