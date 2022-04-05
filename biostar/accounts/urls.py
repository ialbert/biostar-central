from django.contrib import admin
from django.urls import include, path, re_path
from django.conf import settings
from biostar.accounts import views


account_patterns = [
    # Get the reset/ urls
    path('', views.listing, name="accounts_index"),

    path(r'password/reset/', views.password_reset, name='password_reset'),
    path(r'password/reset/done/', views.password_reset_done, name='password_reset_done'),

    path(r'verify/<uidb64>/<token>/', views.email_verify_account, name='email_verify_account'),

    path(r'reset/<uidb64>/<token>/', views.pass_reset_confirm, name='password_reset_confirm'),

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

    # Used for 3rd party logins.
    path("social/", include('allauth.urls')),

]


urlpatterns = [

    path("", include(account_patterns)),

    # Add admin urls.
    path('admin/', admin.site.urls),

]

if settings.PAGEDOWN_IMAGE_UPLOAD_ENABLED:

    urlpatterns += [
        # Pagedown image upload url.
        path('pagedown/image-upload/', views.image_upload_view, name="pagedown-image-upload")
    ]
