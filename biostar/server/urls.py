
from django.conf import settings
from django.urls import path, include
import debug_toolbar
from django.conf.urls.static import static
from biostar.forum.urls import forum_patterns
from biostar.recipes.urls import recipes_patterns
import biostar.accounts.urls as accounts_urls
from biostar.planet.urls import planet_patterns


urlpatterns = [

    # Mount recipes urls on root.
    path(r'', include(recipes_patterns)),

    # Include forum urls
    path(r'forum/', include(forum_patterns)),

    # Include planets urls
    path(r'planet/', include(planet_patterns)),

    # Include the accounts urls
    path(r'accounts/', include(accounts_urls)),

]

if settings.PAGEDOWN_IMAGE_UPLOAD_ENABLED:

    urlpatterns += [
        # Pagedown image upload url.
        path('', include('pagedown.urls'))
    ]


if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT, show_indexes=True)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT, show_indexes=True)
