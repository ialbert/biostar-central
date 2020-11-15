# Import all common settings.
from biostar.settings import *

# Additional apps enabled.
EMAILER_APP = [
    'biostar.emailer.apps.EmailerConfig'
]

INSTALLED_APPS = DEFAULT_APPS + EMAILER_APP

# The url specification.
ROOT_URLCONF = 'biostar.emailer.urls'

# This flag is used flag situation where a data migration is in progress.
# Allows us to turn off certain type of actions (for example sending emails).
DATA_MIGRATION = False

SEND_MAIL = True
