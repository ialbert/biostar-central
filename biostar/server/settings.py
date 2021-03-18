# Used when testing all apps at once

import os

from biostar.settings import *
from biostar.accounts.settings import *
from biostar.emailer.settings import *
from biostar.recipes.settings import *
from biostar.forum.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + PAGEDOWN_APP + PLANET_APPS + ENGINE_APPS + ACCOUNTS_APPS + EMAILER_APP

TASK_RUNNER = 'block'

ROOT_URLCONF = 'biostar.server.urls'
