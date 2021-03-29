from biostar.accounts.settings import *


PLANET_APPS = [
    'biostar.planet.apps.PlanetConfig',

]

# Initialize blog posts upon migration
INIT_PLANET = True

INSTALLED_APPS = PLANET_APPS + ACCOUNTS_APPS

BLOGS_PER_PAGE = 30
PLANET_DIR = os.path.abspath(os.path.join(BASE_DIR, 'export', 'planet'))
