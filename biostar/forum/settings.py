# Inherit from the main settings file.

# Inherit from the accounts settings file.
from biostar.accounts.settings import *

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

SINGLE_FEED_COUNT = 40

SESSION_UPDATE_SECONDS = 40

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

FORUM_APPS = [
    'biostar.forum.apps.ForumConfig',
    'pagedown',
    'debug_toolbar',
]

# Additional middleware.
MIDDLEWARE += [
    'biostar.forum.middleware.forum_middleware',
    'debug_toolbar.middleware.DebugToolbarMiddleware'
]

# Import the default pagedown css first, then our custom CSS sheet
# to avoid having to specify all the default styles
PAGEDOWN_WIDGET_CSS = ('pagedown/demo/browser/demo.css', "lib/pagedown.css",)

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ACCOUNTS_APPS + EMAILER_APP

ROOT_URLCONF = 'biostar.forum.urls'

WSGI_APPLICATION = 'biostar.wsgi.application'

# Time between two accesses from the same IP to qualify as a different view.
POST_VIEW_MINUTES = 7

COUNT_INTERVAL_WEEKS = 10000


try:
    from conf.run.secrets import *
    print(f"Loaded secrets from: conf.run.secrets")
except Exception as exc:
    print(f"Unable to load secrets: {exc}")
