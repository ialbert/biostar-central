# Used when testing all apps at once

from biostar.settings import *
from biostar.accounts.settings import *
from biostar.emailer.settings import *
from biostar.recipes.settings import *
from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + PLANET_APPS + ENGINE_APPS + ACCOUNTS_APPS + EMAILER_APP

ROOT_URLCONF = 'biostar.server.urls'
