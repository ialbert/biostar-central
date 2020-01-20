import os
from biostar.accounts.settings import *


PLANET_APPS = [
    'biostar.planet.apps.PlanetConfig',

]

INSTALLED_APPS = DEFAULT_APPS + PLANET_APPS + ACCOUNTS_APPS + EMAILER_APP

BLOGS_PER_PAGE = 30
PLANET_DIR = os.path.abspath(os.path.join(BASE_DIR, 'biostar', 'planet', 'initial'))


