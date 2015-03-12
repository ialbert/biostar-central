from run.sqlite import *

# add debugging middleware
MIDDLEWARE_CLASSES = list(MIDDLEWARE_CLASSES)
MIDDLEWARE_CLASSES.insert(0, 'debug_toolbar.middleware.DebugToolbarMiddleware')

# We load the debug toolbar as well
INSTALLED_APPS = list(INSTALLED_APPS)

# This needs to be added before the user models.
INSTALLED_APPS.append( "debug_toolbar")

# Enable celery.
CELERY_ENABLED = False

# Must set a broker
BROKER_URL = 'django://'

# Add django transport agent.
# The site must be initialized again.
#INSTALLED_APPS = list(INSTALLED_APPS) +  [ 'kombu.transport.django' ]

# Send a welcome email.
SEND_WELCOME_EMAIL = True

# Send emails via the console.
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'