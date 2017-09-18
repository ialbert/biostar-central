"""website URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
from . import views

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^accounts/login/', views.user_login, name="login"),
    url(r'^signup/$', views.user_signup, name="signup"),
    url(r'^$', views.index, name="index"),
    url(r'^logout/$', views.user_logout, name="logout"),
    url(r'^(?P<id>\d+)/$', views.home, name="home"),
    url(r'^allprojects/$', views.project_list, name='allprojects'),
    url(r'^allprojects/(?P<id>\d+)/$', views.project_detail, name='project'),
    url(r'^allprojects/(?P<id>\d+)/edit/$', views.project_edit, name='editproject'),
    url(r'^allprojects/(?P<id>\d+)/alldata/$', views.data_list, name='alldata'),
    url(r'^allprojects/(?P<id>\d+)/alldata/upload/$', views.data_upload, name='data_upload'),
    url(r'^allprojects/(?P<id>\d+)/allanalysis/$', views.analysis_list, name='allanalysis'),
    url(r'^allprojects/(?P<id>\d+)/alldata/(?P<id2>\d+)/$', views.data_detail, name='data'),
    url(r'^allprojects/(?P<id>\d+)/allanalysis/(?P<id2>\d+)/$', views.analysis_detail, name='analysis'),

]

