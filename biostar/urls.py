from django.conf.urls import url, include
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings


import biostar.accounts.urls as accounts_urls
import biostar.recipes.urls as engine_urls
import biostar.forum.urls as forum_urls

if settings.ENGINE_AS_ROOT:
    engine_url_pattern = r'^'
    forum_url_pattern = r'^forum/'
else:
    engine_url_pattern = r'^projects/'
    forum_url_pattern = r'^'


# Default url patters for the engine.
urlpatterns = [

    # The engine handler.
    url(engine_url_pattern, include(engine_urls)),

    # Forum urls
    url(forum_url_pattern, include(forum_urls)),

    url(r'^accounts/', include(accounts_urls)),

]


if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)
    if settings.INTERNAL_IPS:
        import debug_toolbar
        urlpatterns = [
            url(r'^__debug__/', include(debug_toolbar.urls)),
        ] + urlpatterns


