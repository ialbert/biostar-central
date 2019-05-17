# Import all common settings.
from biostar.accounts.settings import *


# Additional apps enabled.
MESSAGE_APPS = [
    'biostar.message.apps.MessageConfig',
]

INSTALLED_APPS = DEFAULT_APPS + MESSAGE_APPS + ACCOUNTS_APPS

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
                'biostar.message.context.message',
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


LOGIN_REDIRECT_URL = "/inbox/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True


# The url specification.
ROOT_URLCONF = 'biostar.message.urls'

MESSAGES_PER_PAGE = 50

