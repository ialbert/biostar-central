

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^login/', views.user_login, name="login"),
    url(r'^signup/$', views.user_signup, name="signup"),
    url(r'profile/(?P<id>\d+)/$', views.user_profile, name="profile"),
    url(r'^logout/$', views.user_logout, name="logout")
]