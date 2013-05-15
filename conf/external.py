#
# import from the main settings then override some of the values
#
from main.settings import *

# this is the authentication module and function in it
EXTERNAL_AUTHENTICATOR_FUNC = 'dummy_auth'

MIDDLEWARE_CLASSES.extend([
    'main.extauth.ExternalAuthenticator',
])
