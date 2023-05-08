# Inherit from the main settings file.
import os, sys

from biostar.accounts.settings import *

# Inherit from the accounts settings file.
from biostar.planet.settings import *

def join(*args):
    return os.path.abspath(os.path.join(*args))

# Django debug flag.
DEBUG = True

SITE_NAME = 'Biostar Forum'

# Site settings.
POSTS_PER_PAGE = 50
USERS_PER_PAGE = 100
MESSAGES_PER_PAGE = 100
TAGS_PER_PAGE = 50
AWARDS_PER_PAGE = 50

STATS_DIR = os.path.join(BASE_DIR, "export", "stats")


# Enable image upload
PAGEDOWN_IMAGE_UPLOAD_ENABLED = True

# Upload path for pagedown images, relative to media root.
PAGEDOWN_IMAGE_UPLOAD_PATH = "images"

# File path listing tags.
# With one at least being required when making a post.
REQUIRED_TAGS = ''

# Link to display after a post fails to have required tags.
REQUIRED_TAGS_URL = "/"

# How to run tasks in the background.
# Valid options; block, disable, threaded, uwsgi, celery.
TASK_RUNNER = 'threaded'

# Threshold to classify spam
SPAM_THRESHOLD = .5

# Allows post closing.
ALLOW_POST_CLOSING = False

# Classify posts and assign a spam score on creation.
CLASSIFY_SPAM = True

# Log the time for each request
TIME_REQUESTS = True

# Number of results to display in total.
SEARCH_LIMIT = 50

# Initialize the planet app.
INIT_PLANET = False

# Minimum amount of characters to preform searches
SEARCH_CHAR_MIN = 1

# How many posts to index in one job.
BATCH_INDEXING_SIZE = 1000

# Add another context processor to first template.
TEMPLATES[0]['OPTIONS']['context_processors'] += [
    'biostar.forum.context.forum'
]

# Set the number of items in each feed.
VOTE_FEED_COUNT = 7
LOCATION_FEED_COUNT = 7
AWARDS_FEED_COUNT = 7
REPLIES_FEED_COUNT = 15

SIMILAR_FEED_COUNT = 30

SESSION_UPDATE_SECONDS = 10

# Maximum number of awards every SESSION_UPDATE_SECONDS.
MAX_AWARDS = 2

# How many stories to show
HERALD_LIST_COUNT = 100

# Search index name
INDEX_NAME = os.environ.setdefault("INDEX_NAME", "index")
# Relative index directory

INDEX_DIR = os.environ.setdefault("INDEX_DIR", "search")

# Absolute path to index directory in export/
INDEX_DIR = join(MEDIA_ROOT, '..', INDEX_DIR)

# Absolute path to the spam model
SPAM_DATA  = join(BASE_DIR, "export", "spam.data.tar.gz")
SPAM_MODEL = join(BASE_DIR, "export", "spam.model")

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

LOGIN_REDIRECT_URL = "/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

SOCIALACCOUNT_ADAPTER = "biostar.accounts.adapter.SocialAccountAdapter"

FORUM_APPS = [
    'biostar.forum.apps.ForumConfig',
]


VOTE_RATE = '250/h'
EDIT_RATE = '250/h'
SUBS_RATE = '100/h'
DIGEST_RATE = '100/h'

# Additional middleware.
MIDDLEWARE += [
    #'biostar.forum.middleware.ban_ip',
    'biostar.forum.middleware.user_tasks',
    'biostar.forum.middleware.benchmark',
]

# Post types displayed when creating, empty list displays all types.
ALLOWED_POST_TYPES = []


# Import the default pagedown css first, then our custom CSS sheet
# to avoid having to specify all the default styles
PAGEDOWN_WIDGET_CSS = ('pagedown/demo/browser/demo.css',)

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + PAGEDOWN_APP + PLANET_APPS + ACCOUNTS_APPS + EMAILER_APP

# Documentation for docs
FORUM_DOCS = os.path.join(DOCS_ROOT, "forum")

# Add docs to static files directory
STATICFILES_DIRS += [DOCS_ROOT]

ROOT_URLCONF = 'biostar.forum.urls'

WSGI_APPLICATION = 'biostar.wsgi.application'

# Time between two accesses from the same IP to qualify as a different view (seconds)
POST_VIEW_TIMEOUT = 300

# This flag is used flag situation where a data migration is in progress.
# Allows us to turn off certain type of actions (for example sending emails).
DATA_MIGRATION = False

# Default cache
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.dummy.DummyCache',
        #'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake',
    }
}

# Strict rules applied to post tags
STRICT_TAGS = True

DROPDOWN_TAGS = False
TASK_MODULES = ("biostar.forum.tasks", )

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    INSTALLED_APPS.extend([
        'debug_toolbar',
    ])
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')


# Words when present in the content get you banned.
BANNED_WORDS_CONTENT = r"""
\bcialis
\bviagra
\bmoney
\bloan
\bcustomer
\bcash 
"""

# Words, that when present in the title get you banned.
BANNED_WORDS_TITLE  = r"""
\bcash
\bmoney
\bloan
\d{6,}
http 
https
"""


