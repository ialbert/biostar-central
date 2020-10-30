import sys
from biostar.accounts.settings import *
from biostar.emailer.settings import *

# Django debug flag.
DEBUG = True

# Private key used to validate external logins
LOGIN_PRIVATE_KEY = SECRET_KEY

RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_PRIVATE_KEY = ""

SEARCH_CHAR_MIN = 2

SOCIALACCOUNT_EMAIL_VERIFICATION = None
SOCIALACCOUNT_EMAIL_REQUIRED = False
SOCIALACCOUNT_QUERY_EMAIL = True

# Maximum number of projects allowed
MAX_PROJECTS = 20

# Root directory relative to the job path usd to store logs.
JOB_LOGDIR = 'runlog'


PAGEDOWN_IMAGE_UPLOAD_ENABLED = True

# Amount of objects shown per page.
PER_PAGE = 50

# Upload path for pagedown images, relative to media root.
PAGEDOWN_IMAGE_UPLOAD_PATH = "images"

# Stdout filename for each job relative to the log directory.
JOB_STDOUT = os.path.join(JOB_LOGDIR, 'stdout.txt')
JOB_STDERR = os.path.join(JOB_LOGDIR, 'stderr.txt')

# Maximum amount of data allowed
MAX_DATA = 100

# Maximum amount of items per clipboard
MAX_CLIPBOARD = 5

# Name of the clipboard inside of sessions
CLIPBOARD_NAME = "clipboard"

# Maximum amount of total running jobs allowed for non-staff user.
MAX_RUNNING_JOBS = 5

# Maximum amount of cumulative uploaded files a user is allowed, in mega-bytes.
MAX_UPLOAD_SIZE = 10

# Deployment specific parameters.
PROTOCOL = "http"
HTTP_PORT = '8000'
BASE_URL = f"{PROTOCOL}://{SITE_DOMAIN}:{HTTP_PORT}"

FTP_HOST = "localhost"
FTP_PORT = 8021

# Should the site allow signup.
ALLOW_SIGNUP = False

# Allow users to toggle their status moderator
ALLOW_SELF_MODERATE = False

# Maximum size of each file upload in MB
MAX_FILE_SIZE_MB = 300

LOGIN_REDIRECT_URL = "/project/list/"
ACCOUNT_AUTHENTICATED_LOGIN_REDIRECTS = True

ENGINE_APPS = [
    'biostar.recipes.apps.EngineConfig',
    'django.contrib.redirects',
]

PAGEDOWN_APP = ['pagedown.apps.PagedownConfig']

INSTALLED_APPS = DEFAULT_APPS + ENGINE_APPS + ACCOUNTS_APPS + EMAILER_APP + PAGEDOWN_APP

# Additional middleware.
MIDDLEWARE += [
    'biostar.recipes.middleware.recipes_middleware',
]

# The URL configuration.
ROOT_URLCONF = 'biostar.recipes.urls'

# Add another context processor to first template.
TEMPLATES[0]['OPTIONS']['context_processors'] += [
    'biostar.recipes.context.engine'
]

TEMPLATES[0]['OPTIONS']['string_if_invalid'] = "## TEMPLATE ERROR MISSING VARIABLE: %s ##"


# Directory where files are extracted to
EXTRACT_TO = join(BASE_DIR, "export", "extracted")
os.makedirs(EXTRACT_TO, exist_ok=True)

# The location of application specific data.
LOCAL_ROOT = join(BASE_DIR, 'export', 'local')

# The location for the table of contents.
TOC_ROOT = join(MEDIA_ROOT, 'tocs')

IMPORT_ROOT_DIR = join(BASE_DIR, 'export', 'local')

# Configure language detection
LANGUAGE_DETECTION = ['en']

# Ensure that the table of directory exists.
os.makedirs(TOC_ROOT, exist_ok=True)

# Sendfile settings go here.
SENDFILE_ROOT = MEDIA_ROOT
SENDFILE_URL = '/protected/'

COUNT_INTERVAL_WEEKS = 10000

SESSION_ENGINE = "django.contrib.sessions.backends.db"

# SENDFILE_BACKEND = "sendfile.backends.nginx"

SENDFILE_BACKEND = "sendfile.backends.development"

STATICFILES_FINDERS += [
    'compressor.finders.CompressorFinder',
]

# Tries to load up secret settings from a predetermined module

try:
    from conf.run.secrets import *
    #print(f"Loaded secrets from: conf.run.secrets")
except Exception as exc:
    #print(f"Secrets module not imported: {exc}", file=sys.stderr)
    pass


# Enable debug toolbar specific functions
if DEBUG and DEBUG_TOOLBAR and 'debug_toolbar' not in INSTALLED_APPS:
    INSTALLED_APPS.extend(['debug_toolbar'])
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')
