
from django.conf import settings
from django.urls import path, include

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
