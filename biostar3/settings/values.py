#
# Site specic behaviors. These values are loaded first even before the base settings module.
#
from collections import OrderedDict

# Should the site send a welcome email.
SEND_WELCOME_EMAIL = True


# Google ReCaptcha No-Captcha settings
# When set the captcha forms will be active.
RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_SECRET_KEY = ""

# Enable rate limiting.
RATELIMIT_ENABLE = True

# How many maximum signup accesses per minute.
# See django-ratelimit for rates and keys.
# Does not include signups via social authentication.
SIGNUP_RATELIMIT = "3/m"

# How many recent votes to show.
RECENT_VOTE_COUNT = 10

# How many recent users to show.
RECENT_USER_COUNT = 10

# How many posts per page.
POST_PAGINATE_BY = 10

# The CSS classes associated with the Django messages framework.
MESSAGE_TAGS = {
    10: 'alert-info', 20: 'alert-info', 25: 'alert-success', 30: 'alert-warning', 40: 'alert-danger',
}

# Connects a word in post sort to a queryset sort attribute of the data model.
POST_SORT_MAP = OrderedDict([
    ("update", "-lastedit_date"),
    ("views", "-view_count"),
    ("followers", "-subs_count"),
    ("answers", "-reply_count"),
    ("bookmarks", "-book_count"),
    ("votes", "-vote_count"),
    ("rank", "-rank"),
    ("creation", "-creation_date"),
])

# These are the fields rendered in the post sort order drop down.
POST_SORT_FIELDS = POST_SORT_MAP.keys()
POST_SORT_DEFAULT = POST_SORT_FIELDS[0]

# The messages show when the sort is not valid.
POST_SORT_INVALID_MSG = "Invalid sort parameter in URL."

# Messaging related settings.
LOCAL_MESSAGE, EMAIL_MESSAGE, NO_MESSAGES, SMART_MESSAGES, ALL_MESSAGES = range(5)

MESSAGING_MAP = OrderedDict([
    (SMART_MESSAGES, "smart mode",),
    (LOCAL_MESSAGE, "local messages",),
    (EMAIL_MESSAGE, "emails",),
    (ALL_MESSAGES, "email for every new thread (mailing list mode)",),
])

MESSAGE_DEFAULT = SMART_MESSAGES

#
# Django allauth configuration
#
SOCIALACCOUNT_ADAPTER = 'biostar3.middleware.AutoSignupAdapter'

ACCOUNT_AUTHENTICATION_METHOD = "email"
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_UNIQUE_EMAIL = True
ACCOUNT_USERNAME_REQUIRED = False
ACCOUNT_EMAIL_VERIFICATION = "optional"
ACCOUNT_EMAIL_SUBJECT_PREFIX = "[biostar]"
ACCOUNT_USER_DISPLAY = lambda user: user.name
ACCOUNT_PASSWORD_MIN_LENGHT = 6
ACCOUNT_USER_MODEL_USERNAME_FIELD = None
ACCOUNT_USER_MODEL_EMAIL_FIELD = "email"
ACCOUNT_DEFAULT_HTTP_PROTOCOL = "http"
ACCOUNT_CONFIRM_EMAIL_ON_GET = True
LOGIN_REDIRECT_URL = "/"
ACCOUNT_SESSION_COOKIE_AGE = 3600 * 24 * 60