# Import all common settings.
from biostar.accounts.settings import *


# Additional apps enabled.
INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    'mailer',
    'compressor',
    'taggit',
    'snowpenguin.django.recaptcha2',
    'rest_framework',

    # The order of apps matters in the template loading
    'biostar.message.apps.MessageConfig',
    'biostar.accounts.apps.AccountsConfig',

    # Allauth templates come last.
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'allauth.socialaccount.providers.github',
]


# The url specification.
ROOT_URLCONF = 'biostar.message.urls'

MESSAGES_PER_PAGE = 50

