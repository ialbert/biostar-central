from biostar.settings import *

import os

# Apply the logger settings.
from biostar.logconf import LOGGING

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# Django debug flag.
DEBUG = True

# Should the site allow signup.
ALLOW_SIGNUP = False

# Private key used to validate external logins
LOGIN_PRIVATE_KEY = "private-key"

# In MB
MAX_UPLOAD_SIZE = 10

# The password for admin users. Must be changed in production.
DEFAULT_ADMIN_PASSWORD = "admin@localhost"

RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_PRIVATE_KEY = ""

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True


# Should the site allow signup.
ALLOW_SIGNUP = False

LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

INSTALLED_APPS += [

    # Accounts configuration.
    'biostar.accounts.apps.AccountsConfig',

    # Allauth templates come last.
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'allauth.socialaccount.providers.github',
]

# Template specific settings.
TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'APP_DIRS': True,
        'OPTIONS': {
            'string_if_invalid': "**MISSING**",
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.template.context_processors.media',
                'django.contrib.messages.context_processors.messages',
                'biostar.accounts.context.accounts',
            ],
            # 'loaders': [
            #     ('django.template.loaders.cached.Loader',
            #         'django.template.loaders.filesystem.Loader',
            #         'django.template.loaders.app_directories.Loader',
            #     )
            # ]
        },
    },
]


ROOT_URLCONF = 'biostar.accounts.urls'

# List of social login clients tuples.
# ( name, client_id, secret )

# Default clients redirect to localhost.
# Default clients may not be operational. See the
# django allauth documentation on how to set them up.

#
# Callback example settings:
#
# http://localhost:8000/accounts/social/google/login/callback/
# http://localhost:8000/accounts/social/github/login/callback/
#
SOCIAL_CLIENTS = [

    ("Google", "547073349197-ri3eku9fdpi1ble7eoc4amlsrh7m2oiv.apps.googleusercontent.com", "DR1-zMqOLTqRGvhSDb5rQBMg"),
    ("GitHub", "d8493ce8967ea5abbd73", "04ff5043ebce2317d3c26bfe90fbc5e67fa38d05")

]
