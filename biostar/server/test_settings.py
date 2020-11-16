from .settings import *

# Do not multi-thread tests.
MULTI_THREAD = False

# Skip hitting the spam indexe when creating test posts
CLASSIFY_SPAM = False


# Turn the emailing tasks off for tests
SEND_MAIL = False