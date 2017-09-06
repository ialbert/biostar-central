from django.conf.urls import url
from . import views
from .forms import CustomLoginForm

urlpatterns = [
    url(r'^$', views.index, name="home"),
    url(r'^login/$', views.auth_views.login, 
    {"template_name":"login.html", "authentication_form": CustomLoginForm}, name="login"),
    url(r'^signup/$', views.signup, name="signup"),
    url(r'^logout/$', views.auth_views.logout, name="logout"),
]
