# Import all common settings.
from biostar.settings import *

# Additional apps enabled.
EMAILER_APP = [
    'biostar.emailer.apps.EmailerConfig'
]

INSTALLED_APPS = DEFAULT_APPS + EMAILER_APP

# The url specification.
ROOT_URLCONF = 'biostar.emailer.urls'


