from django.contrib import admin
from django.urls import include, path, re_path
from django.conf import settings
from biostar.accounts import views


account_patterns = [
    # Get the reset/ urls
    path('', views.listing, name="accounts_index"),
    path('admin/', admin.site.urls, name='django_admin'),

    path(r'password/reset/', views.password_reset, name='password_reset'),
    path(r'password/reset/done/', views.password_reset_done, name='password_reset_done'),

    re_path(r'verify/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/',
            views.email_verify_account, name='email_verify_account'),

    re_path(r'reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/',
            views.pass_reset_confirm, name='password_reset_confirm'),

    path(r'verify/', views.send_email_verify, name="send_email_verify"),

    path(r'reset/done/', views.password_reset_complete, name='password_reset_complete'),
    path(r'moderate/<int:uid>/', views.user_moderate, name="user_moderate"),
    path(r'login/', views.user_login, name="login"),
    path(r'signup/', views.user_signup, name="signup"),
    path(r'profile/<str:uid>/', views.user_profile, name="user_profile"),

    path(r'edit/profile/', views.edit_profile, name='edit_profile'),
    path(r'toggle/notify/', views.toggle_notify, name='toggle_notify'),
    path(r'logout/', views.user_logout, name="logout"),

    path(r'debug/user/', views.debug_user, name="debug_user"),

    # Message urls
    path(r'inbox/', views.message_list, name='inbox'),

    # External url login
    path(r'external/', views.external_login, name="external"),

    # Used for 3rd party logins.
    path("social/", include('allauth.urls')),

]


urlpatterns = [

    path("", include(account_patterns)),

]
