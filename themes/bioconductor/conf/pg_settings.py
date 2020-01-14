from themes.bioconductor.settings import *

INSTALLED_APPS = DEFAULT_APPS + FORUM_APPS + ACCOUNTS_APPS + EMAILER_APP

DEBUG = True

# Show debug toolbar
DEBUG_TOOLBAR = DEBUG

# Enable debug toolbar specific functions
if DEBUG_TOOLBAR:
    INSTALLED_APPS.extend([
        'debug_toolbar',
    ])
    MIDDLEWARE.append('debug_toolbar.middleware.DebugToolbarMiddleware')

# Add static serving capability
if DEBUG is False:
    print("Whitenoise static serve enabled (pip install whitenoise)")
    MIDDLEWARE.append('whitenoise.middleware.WhiteNoiseMiddleware')
