from biostar.settings import *

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
    'pagedown',
    'taggit',
    'debug_toolbar',
    'snowpenguin.django.recaptcha2',
    'rest_framework',
    # The order of apps matters in the template loading

    'biostar.forum.apps.ForumConfig',
    'biostar.emailer.apps.EmailerConfig',
    'biostar.accounts.apps.AccountsConfig',
    'biostar.message.apps.MessageConfig',

    # Allauth templates come last.
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'allauth.socialaccount.providers.github',
]
