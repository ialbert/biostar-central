
from django.conf import settings
from django.conf.urls import url, include

import biostar.forum.urls as forum_urls
import biostar.recipes.urls as engine_urls
import biostar.accounts.urls as accounts_urls


urlpatterns = [

    url(r'^', include(forum_urls)),

    # Include messages urls
    url(r'^recipes/', include(engine_urls)),

    # Include the accounts urls
    # Already included with above urls.
    #url(r'^accounts/', include(accounts_urls)),

]

    #print(p.url_patterns)

#print(urlpatterns)