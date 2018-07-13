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

    # The user account handler.
    ACCOUNTS,

    # The django generated admin site.
    ADMIN

]

# Have the option to load the engine and forum together
if settings.ENABLE_FORUM:
    urlpatterns += [url(r'^forum/', include(forum_urls))]

# Urls mounted when forum is enabled by itself
if settings.ONLY_ENABLE_FORUM:

    urlpatterns = [

                # The forum handler.
                url(r'^', include(forum_urls)),

                # The user account handler.
                ACCOUNTS,

                # The django generated admin site.
                ADMIN,
                ]


# Add the message urls
urlpatterns += forum_urls.msg_urls

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)

