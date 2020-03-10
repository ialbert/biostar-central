from biostar.server.settings import *

# Turn multi threading off during testing.
MULTI_THREAD = False

# reCaptcha left alone during testing
RECAPTCHA_PUBLIC_KEY = ''
RECAPTCHA_PRIVATE_KEY = ''

# Don't hammer away at the rss feeds while testing.
INIT_PLANET = False
