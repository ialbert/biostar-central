# Import all common settings.
from biostar.settings import *

# Additional apps enabled.
INSTALLED_APPS += [
    'biostar.emailer.apps.EmailerConfig'
]

# The url specification.
ROOT_URLCONF = 'biostar.emailer.urls'


