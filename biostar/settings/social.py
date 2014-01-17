# Settings for social authentication.

__author__ = 'ialbert'

import os
from django.core.exceptions import ImproperlyConfigured

def get_env(name):
    """Get the environment variable or return exception"""
    try:
        return os.environ[name]
    except KeyError:
        msg = "*** Required environment variable %s not set." % name
        raise ImproperlyConfigured(msg)


AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',
    "allauth.account.auth_backends.AuthenticationBackend",
)

SOCIALACCOUNT_PROVIDERS = {
    #'facebook': {
    #    'SCOPE': ['email'],
    #    'AUTH_PARAMS': {'auth_type': 'reauthenticate'},
    #    'METHOD': 'oauth2',
    #    'LOCALE_FUNC': lambda x: 'en_US',
    #},
    'google': {
        'SCOPE': ['email', 'https://www.googleapis.com/auth/userinfo.profile'],
        'AUTH_PARAMS': {'access_type': 'online'},
        'PROVIDER_KEY' : get_env("GOOGLE_PROVIDER_KEY"),
        'PROVIDER_SECRET_KEY': get_env("GOOGLE_PROVIDER_SECRET_KEY"),
    },
}

ACCOUNT_AUTHENTICATION_METHOD = "email"
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_UNIQUE_EMAIL = True
ACCOUNT_USERNAME_REQUIRED = False
ACCOUNT_EMAIL_VERIFICATION = "optional"
ACCOUNT_EMAIL_SUBJECT_PREFIX = "[biostar] "
ACCOUNT_PASSWORD_MIN_LENGHT = 6
ACCOUNT_USER_MODEL_USERNAME_FIELD="email"

ACCOUNT_DEFAULT_HTTP_PROTOCOL = "http"
