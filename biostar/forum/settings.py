# Inherit from the main settings file.
import os, sys
from biostar.accounts.settings import *

# Inherit from the accounts settings file.
from biostar.planet.settings import *

# Django debug flag.
DEBUG = True

SITE_NAME = 'Biostar Forum'

# Site settings.
POSTS_PER_PAGE = 50
USERS_PER_PAGE = 100
MESSAGES_PER_PAGE = 100
TAGS_PER_PAGE = 50

STATS_DIR = os.path.join(BASE_DIR, "export", "stats")

REQUIRED_TAGS = ""

# Time period to cache Ips for banning.
TIME_PERIOD = 24 * 3600

# How many visit within that time period.
MAX_VISITS = 50

# Whitelist of Ip addresses.
IP_WHITELIST = [

]


PAGEDOWN_IMAGE_UPLOAD_ENABLED = True

# Upload path for pagedown images, relative to media root.
PAGEDOWN_IMAGE_UPLOAD_PATH = "images"

REQUIRED_TAGS_URL = "/"

BANNED_IPS = os.path.join(BASE_DIR, 'export', 'logs', 'banned.txt')


# Posts cut off when applying a filtre
CUTOFF = 1000

# The gravatar image used for users, applied to all users.
GRAVATAR_ICON = ''

SPAM_THRESHOLD = .5

# Spam index used to classify new posts as spam or ham.
SPAM_INDEX_NAME = os.getenv("SPAM_INDEX_NAME", "spam")

SPAM_INDEX_DIR = 'spammers'

# Absolute path to spam index directory in export/
SPAM_INDEX_DIR = os.path.abspath(os.path.join(MEDIA_ROOT, '..', SPAM_INDEX_DIR))

# Classify posts and assign a spam score on creation.
CLASSIFY_SPAM = True

ENABLE_DIGESTS = False

# Disable all asynchronous tasks
DISABLE_TASKS = False

# Log the time for each request
TIME_REQUESTS = True

# Indexing interval in seconds.
INDEX_SECS_INTERVAL = 10

# Number of results to display in total.
SEARCH_LIMIT = 20

INIT_PLANET = False

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

PAGEDOWN_APP = ['pagedown.apps.PagedownConfig']

FORUM_APPS = [

    'biostar.forum.apps.ForumConfig',
]

# Additional middleware.
MIDDLEWARE += [
    'biostar.forum.middleware.ban_ip',
    'biostar.forum.middleware.user_tasks',
    'biostar.forum.middleware.benchmark',
]

# Remap the post type display to a more human friendly one.
REMAP_TYPE_DISPLAY = False

# Post types displayed when creating, empty list displays all types.
ALLOWED_POST_TYPES = []


# Import the default pagedown css first, then our custom CSS sheet
# to avoid having to specify all the default styles
PAGEDOWN_WIDGET_CSS = ('pagedown/demo/browser/demo.css',)

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + PAGEDOWN_APP + PLANET_APPS + ACCOUNTS_APPS + EMAILER_APP

FORUM_DOCS = os.path.join(DOCS_ROOT, "forum")

# Add docs to static files directory
STATICFILES_DIRS += [DOCS_ROOT]


# Directory for the planets app.
#PLANET_DIR = ''

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
    #print(f"Loaded secrets from: conf.run.secrets")
except Exception as exc:
    print(f"Secrets module not imported: {exc}", file=sys.stderr)
    pass

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    INSTALLED_APPS.extend([
        'debug_toolbar',
    ])
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')
