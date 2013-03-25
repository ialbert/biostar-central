"""
How to change the settings for a different site
"""

# import from the main settings then override some of them
from main.settings import *

DEBUG = 1

TEMPLATE_DEBUG = DEBUG

# admin site may fail if this setting is active
#TEMPLATE_STRING_IF_INVALID = "***"

# set the domain administrators
ADMINS = (
    ('Istvan Albert', 'istvan.albert@gmail.com'),
)

# set the site url
SITE_DOMAIN = 'localhost:8080'

# set the secret key for the site
SECRET_KEY = 'murjkj468712u7u2888271209239929u7u2888271209239929u7u28882'

SESSION_UPDATE_TIME = 600

SELENIUM_TEST_LOGIN_TOKEN = "abc"

EMAIL_HOST = 'smtp.psu.edu'
EMAIL_HOST_USER = ''
EMAIL_HOST_PASSWORD = ''
DEFAULT_FROM_EMAIL = 'iua1@psu.edu'
SERVER_EMAIL = 'iua1@psu.edu'

# add a custom html page
TEMPLATE_ROWS['job'] = "row.job-sponsored.html"
LIVE_DIR = path(HOME_DIR, '..', 'live')
TEMPLATE_DIRS.append(LIVE_DIR)
    