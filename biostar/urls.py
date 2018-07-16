from django.conf.urls import url, include
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings

import biostar.accounts.urls as accounts_urls
import biostar.engine.urls as engine_urls
import biostar.forum.urls as forum_urls


ACCOUNTS = url(r'^accounts/', include(accounts_urls))
ADMIN = url(r'^django/admin/', admin.site.urls, name='django_admin')


# Default url patters for the engine.
urlpatterns = [

    # The engine handler.
    url(r'^', include(engine_urls)),

    # The django generated admin site.
    ADMIN,

    ACCOUNTS,

    # Add the message urls at the end
    *forum_urls.msg_urls

]

# Have the option to load the engine and forum together
if settings.ENABLE_FORUM:
    urlpatterns += [url(r'^forum/', include(forum_urls))]

# Urls mounted when forum is enabled by itself
if settings.ONLY_FORUM_URLS:

    # Replace the engine handler with the forums
    urlpatterns[0] = url(r'^', include(forum_urls))


if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)

