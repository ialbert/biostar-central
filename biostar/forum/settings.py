# Inherit from the main settings file.

# Inherit from the accounts settings file.
from biostar.accounts.settings import *
from biostar.message.settings import *

# Django debug flag.
DEBUG = True

# Override compression if needed.
# COMPRESS_ENABLED = True

POSTS_PER_PAGE = 40
USERS_PER_PAGE = 100
MESSAGES_PER_PAGE = 100
TAGS_PER_PAGE = 50

VOTE_FEED_COUNT = 10
LOCATION_FEED_COUNT = 5
AWARDS_FEED_COUNT = 10
REPLIES_FEED_COUNT = 15

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

FORUM_APPS = [
    'biostar.forum.apps.ForumConfig',
]

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + MESSAGE_APPS + ACCOUNTS_APPS

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
                'biostar.forum.context.forum',
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

ROOT_URLCONF = 'biostar.forum.urls'

WSGI_APPLICATION = 'biostar.wsgi.application'


# Time between two accesses from the same IP to qualify as a different view.
POST_VIEW_MINUTES = 7

COUNT_INTERVAL_WEEKS = 10000

