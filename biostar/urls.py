from django.conf.urls import url, include
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings

import biostar.accounts.urls as accounts_urls
import biostar.engine.urls as engine_urls

urlpatterns = [

    # The engine handler.
    url(r'^', include(engine_urls)),

    # The user account handler.
    url(r'^accounts/', include(accounts_urls)),

    # The django generated admin site.
    url(r'^django/admin/', admin.site.urls, name='django_admin'),

]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)

