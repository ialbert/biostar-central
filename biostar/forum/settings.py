# Inherit from the main settings file.

# Inherit from the accounts settings file.
from biostar.accounts.settings import *

# Django debug flag.
DEBUG = True

SITE_NAME = 'Biostar Forum'

# Show debug toolbar
DEBUG_TOOLBAR = False

# Override compression if needed.
# COMPRESS_ENABLED = True

POSTS_PER_PAGE = 40
USERS_PER_PAGE = 100
MESSAGES_PER_PAGE = 100
TAGS_PER_PAGE = 50


# User above this score do not get a reCAPTCHA
RECAPTCHA_TRUSTED_USER_SCORE = 15

# Full file path to tags.
# Each line is a list of comma separated tags.
TAGS_OPTIONS_FILE = ''


# Add list of contributors to the post list.
ADD_THREAD_USERS = True

# The gravatar image used for users, applied to all users.
GRAVATAR_ICON = ''

ENABLE_DIGESTS = False

# Log the time for each request
TIME_REQUESTS = True

# Indexing interval in seconds.
INDEX_SECS_INTERVAL = 10
# Number of results to display in total.
SEARCH_LIMIT = 20

# Minimum amount of characters to preform searches
SEARCH_CHAR_MIN = 1

# Number of results to display per page.
SEARCH_RESULTS_PER_PAGE = 50

BATCH_INDEXING_SIZE = 1000

# Add another context processor to first template.
TEMPLATES[0]['OPTIONS']['context_processors'] += [
    'biostar.forum.context.forum'
]

VOTE_FEED_COUNT = 10
LOCATION_FEED_COUNT = 5
AWARDS_FEED_COUNT = 10
REPLIES_FEED_COUNT = 15

SIMILAR_FEED_COUNT = 30

SESSION_UPDATE_SECONDS = 40

# Search index name
INDEX_NAME = os.environ.setdefault("INDEX_NAME", "index")
# Relative index directory
INDEX_DIR = os.environ.setdefault("INDEX_DIR", "search")
# Absolute path to index directory in export/
INDEX_DIR = os.path.abspath(os.path.join(MEDIA_ROOT, '..', INDEX_DIR))

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

FORUM_APPS = [
    'biostar.forum.apps.ForumConfig',
    'pagedown',
    'django_elasticsearch_dsl',
    'django_elasticsearch_dsl_drf'

]

# Additional middleware.
MIDDLEWARE += [
    'biostar.forum.middleware.user_tasks',
    'biostar.forum.middleware.benchmark',
]

# Remap the post type display to a more human friendly one.
REMAP_TYPE_DISPLAY = False
# Post types displayed when creating, empty list displays all types.
ALLOWED_POST_TYPES = []

ELASTICSEARCH_DSL = {
    'default': {
        'hosts': 'localhost:9200'
    },
}

USE_ELASTIC_SEARCH = False

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    FORUM_APPS.append('debug_toolbar')
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')


# Import the default pagedown css first, then our custom CSS sheet
# to avoid having to specify all the default styles
PAGEDOWN_WIDGET_CSS = ('pagedown/demo/browser/demo.css', "lib/pagedown.css",)

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ACCOUNTS_APPS + EMAILER_APP

ROOT_URLCONF = 'biostar.forum.urls'

WSGI_APPLICATION = 'biostar.wsgi.application'

# Time between two accesses from the same IP to qualify as a different view.
POST_VIEW_MINUTES = 7

COUNT_INTERVAL_WEEKS = 10000

# This flag is used flag situation where a data migration is in progress.
# Allows us to turn off certain type of actions (for example sending emails).
DATA_MIGRATION = False

# Tries to load up secret settings from a predetermined module
# This is for convenience only!
try:
    from conf.run.secrets import *
    print(f"Loaded secrets from: conf.run.secrets")
except Exception as exc:
    print(f"Secrets module not imported: {exc}")

