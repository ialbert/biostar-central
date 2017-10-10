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
from django.conf.urls.static import static
from django.conf import settings
from . import views, user_views

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^accounts/login/', user_views.user_login, name="login"),
    url(r'^signup/$', user_views.user_signup, name="signup"),
    url(r'^$', views.index, name="index"),
    url(r'^logout/$', user_views.user_logout, name="logout"),
    url(r'^project/create/$', views.project_create, name='project_create'),
    url(r'^project/list/$', views.project_list, name='project_list'),
    url(r'^analysis/list/(?P<id>\d+)/$', views.analysis_list, name='analysis_list'),
    url(r'^project/view/(?P<id>\d+)/$', views.project_view, name='project_view'),
    url(r'^project/edit/(?P<id>\d+)/$', views.project_edit, name='project_edit'),
    url(r'^data/list/(?P<id>\d+)/$', views.data_list, name='data_list'),
    url(r'^data/view/(?P<id>\d+)/$', views.data_view, name='data_view'),
    url(r'^data/edit/(?P<id>\d+)/$', views.data_edit, name='data_edit'),
    url(r'^data/create/(?P<id>\d+)/$', views.data_create, name='data_create'),
    url(r'^analysis/view/(?P<id>\d+)/(?P<id2>\d+)/$', views.analysis_view, name='analysis_view'),
    url(r'^analysis/run/(?P<id>\d+)/(?P<id2>\d+)/$', views.analysis_run, name='analysis_run'),
    url(r'^analysis/edit/(?P<id>\d+)/(?P<id2>\d+)/$', views.analysis_edit, name='analysis_edit'),
    url(r'^jobs/list/(?P<id>\d+)/$', views.jobs_list, name='jobs_list'),
    url(r'^job/view/(?P<id>\d+)/$', views.job_view, name='job_view')
]


urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)