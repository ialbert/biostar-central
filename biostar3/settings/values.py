#
# Site specic behaviors.
# These values are loaded first even before the base settings module.
#
# These are located in settings to allow overriding them.
#
from .celeryconfig import *
from collections import OrderedDict


# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True
TEMPLATE_DEBUG = True

# What domain to set the cookies for.
SESSION_COOKIE_DOMAIN = ".lvh.me"

# Site administrators. Make sure to override this.
ADMINS = (
    ("Biostar Community", "1@localhost.com"),
)

# SECURITY WARNING: make this private.
# By default same as the email.
SECRET_KEY = "1@localhost.com"

# The password that needs to be sent to update flat pages.
FLATPAGE_SECRET_KEY = "foo"

MANAGERS = ADMINS

# Site setup.
# lvh.me resolves to localhost and contains dotted domain name.
SITE_ID = 1
SITE_NAME = "Site Name"
SITE_DOMAIN = "www.lvh.me:8080"
SITE_SCHEME = "http"

# Subdomain relative to site domain.
SITE_PREPEND_SUBDOMAIN = False

# Which subdomains lead to the main site.
DEFAULT_SUBDOMAINS = { 'www', '127' }

# Default group name.
DEFAULT_GROUP_NAME = "Biostar"

# Default group domain. Don't change it after the site is deployed.
DEFAULT_GROUP_DOMAIN = "www"

# This must be set correctly in production.
ALLOWED_HOSTS = [".lvh.me"]

# Should the site send a welcome email.
SEND_WELCOME_EMAIL = True

# The email pattern to reply to.
EMAIL_ADDRESS_PATTERN = "{}+{}+{}@lvh.me"

# The email secret key required to process the incoming emails.
EMAIL_HANDLER_SECRET_KEY = "foo"

# Should we remove the quoted text.
EMAIL_REPLY_REMOVE_QUOTED_TEXT = True

# Should the site allow interacting via email.
ALLOW_EMAIL_REPLY = True



# Google ReCaptcha No-Captcha settings
# When set the captcha forms will be active.
RECAPTCHA_PUBLIC_KEY = ""
RECAPTCHA_SECRET_KEY = ""

# Enable rate limiting.
RATELIMIT_ENABLE = True

# Which email backend to use.
# Backends may require extra parameters.
EMAIL_BACKEND = 'django.core.mail.backends.locmem.EmailBackend'

# Maximal post size in characters
MAX_POST_SIZE = 150000

# How many maximum signup accesses per minute.
# See django-ratelimit for rates and keys.
# Does not include signups via social authentication.
SIGNUP_RATELIMIT = "30/h"

# How many recent votes to show.
RECENT_VOTE_COUNT = 7

# How many recent users to show.
RECENT_USER_COUNT = 7

# How many recent awards to show.
RECENT_AWARD_COUNT = 7

# How frequently to check for updates.
SESSION_UPDATE_SECONDS = 600

# How many posts per page.
POSTS_PER_PAGE = 10

# How many minutes until a post view from an IP is counted again.
POST_VIEW_INTERVAL = 5


# 10MB -> 10 * 1024 * 1024
MAX_UPLOAD_SIZE = 15 * 1024 * 1024



# Html sanitization. Whitelisting the allowed html content for bleach.sanitize
# Moderators will be allowed to use ALLOWED + TRUSTED settings.
ALLOWED_TAGS = "p div br code pre h1 h2 h3 h4 hr span s sub sup b i img strong \
    strike em underline super table thead tr th td tbody".split()
TRUSTED_TAGS = "embed".split()

ALLOWED_STYLES = 'color font-weight background-color width height'.split()
TRUSTED_STYLES = 'div'.split()

ALLOWED_ATTRIBUTES = {
    '*': ['class', 'style'],
    'a': ['href', 'rel'],
    'img': ['src', 'alt', 'width', 'height'],
    'table': ['border', 'cellpadding', 'cellspacing'],

}
TRUSTED_ATTRIBUTES = {

}

# Secret keys that allows other sites to submit content
FEDERATION_SECRET_KEYS = {
    "foo": ("http://www.foo.com", "Foo Site", "foo-secret"),
    "bar": ("http://www.bar.com", "Bar Site", "bar-secret"),
}

# Which social account provider to trust with sending the correct email.
TRUSTED_SOCIALACCOUNT_PROVIDERS = {
    'google',
    'github',
    'persona',
}


def GET_SUBDOMAIN(request):
    """
    Used to extract the subdomain.
    The front end server must correctly forward the HTTP_HOST header.
    """
    domain = request.META.get('HTTP_HOST', 'www')
    subdomain = domain.split('.')[0]
    return subdomain.lower()

# The value, url pairs that show up
# in the shortcuts widget. The name may
# contain html.
DEFAULT_SHORTCUTS = [
    # Item name, Item URL
    ("Home", "/"),
    ("Unanswered", "/p/unanswered/"),
    ("Jobs", "/t/jobs/"),
    ("All tags", "/t/"),
    ("<b>Your Account</b>", "/site/me/"),
    ("<b>Edit Profile</b>", "/site/edit/my/profile/"),
]

# How many groups can a regular user create.
GROUP_COUNT_PER_USER = 2

# Minimum reputation to create a group.
GROUP_MIN_SCORE = 10000

# Sort values for userlist.
USER_SORT_BY_VISIT, USER_SORT_BY_REP = "visit", "reputation"
USER_SORT_BY_JOIN, USER_SORT_BY_ACTIVITY = "join", "activity"

# Interface representation for user list sort values.
USER_SORT_CHOICES = [
    (USER_SORT_BY_VISIT, "sort by recent visit"),
    (USER_SORT_BY_REP, "sort reputation"),
    (USER_SORT_BY_JOIN, "sort date joined"),
    (USER_SORT_BY_ACTIVITY, "sort activity level"),
]

# Connecting a value to a order by clause on the database.
USER_SORT_ORDER = {
    USER_SORT_BY_VISIT: "-profile__last_login",
    USER_SORT_BY_REP: "-score",
    USER_SORT_BY_JOIN: "profile__date_joined",
    USER_SORT_BY_ACTIVITY: "-activity",
}

# Quick lookup for valid values.
USER_SORT_MAP = dict(USER_SORT_CHOICES)

# The default user sort order.
USER_SORT_DEFAULT = USER_SORT_BY_VISIT

# Error message when passing an incorrect sort parameter.
USER_SORT_INVALID_MSG = "Invalid user sort parameter received"

# Used to generate the sort dropdown.
SORT_BY_UPDATE, SORT_BY_VIEWS = "update", "views"
SORT_BY_SUBS, SORT_BY_ANSWERS, SORT_BY_BOOKMARKS = "followers", "replies", "bookmarks"
SORT_BY_VOTES, SORT_BY_RANK, SORT_BY_CREATION = "votes", "rank", "creation"

# User interface elements for post sort choices.
POST_SORT_CHOICES = [
    (SORT_BY_UPDATE, "sort by activity",),
    (SORT_BY_VIEWS, "sort by views"),
    (SORT_BY_SUBS, "sort by followers"),
    (SORT_BY_ANSWERS, "sort by answers"),
    (SORT_BY_BOOKMARKS, "sort by bookmarks"),
    (SORT_BY_VOTES, "sort by votes"),
    (SORT_BY_RANK, "sort by rank"),
    (SORT_BY_CREATION, "sort by creation"),
]

# Default sort order.
POST_SORT_DEFAULT = SORT_BY_UPDATE

# Quicker lookup of value to label.
POST_SORT_MAP = dict(POST_SORT_CHOICES)

# Connects a sort value to an order_by attribute in the database.
POST_SORT_ORDER = {
    SORT_BY_UPDATE: "-lastedit_date",
    SORT_BY_VIEWS: "-view_count",
    SORT_BY_SUBS: "-subs_count",
    SORT_BY_ANSWERS: "-reply_count",
    SORT_BY_BOOKMARKS: "-book_count",
    SORT_BY_VOTES: "-vote_count",
    SORT_BY_RANK: "-rank",
    SORT_BY_CREATION: "-creation_date"
}

# The messages show when the sort is not valid.
POST_SORT_INVALID_MSG = "Invalid sort parameter in URL."

# Messaging related settings.
LOCAL_TRACKER, EMAIL_TRACKER, LOCAL_MESSAGES, DEFAULT_MESSAGES, MAILING_LIST_MESSAGES, LEAVE_GROUP = range(6)

# The mapping from a messaging type to a readable word.
MESSAGE_CHOICES = [
    (DEFAULT_MESSAGES, "Smart mode"),
    (LOCAL_MESSAGES, "Local messages"),
    (LOCAL_TRACKER, "Local tracker"),
    (EMAIL_TRACKER, "Email tracker"),
    (MAILING_LIST_MESSAGES, "Mailing list"),
    (LEAVE_GROUP, "Leave Group"),
]

EMAIL_MESSAGE_TYPES = {
    EMAIL_TRACKER, MAILING_LIST_MESSAGES,
}

# Default messaging value for a new user.
MESSAGE_DEFAULT = DEFAULT_MESSAGES

# TODO refactoring here
SUBSCRIPTION_CHOICES = MESSAGE_CHOICES
SUBSCRIPTION_DEFAULT = DEFAULT_MESSAGES

# Connects a word to a number of days. 0 indicates no limit.
ALL_TIME, THIS_DAY, THIS_WEEK, THIS_MONTH, THIS_YEAR = 0, 1, 7, 36, 365

TIME_LIMIT_CHOICES = [
    (ALL_TIME, "all time"),
    (THIS_DAY, "today"),
    (THIS_WEEK, "this week"),
    (THIS_MONTH, "this month"),
    (THIS_YEAR, "this year"),
]

# The default time limit value.
TIME_LIMIT_DEFAULT = ALL_TIME

TIME_LIMIT_MAP = dict(TIME_LIMIT_CHOICES)

# Error message for invalid time limit.
TIME_LIMIT_INVALID_MSG = "Invalid time limit received"

# Digest choices.
NO_DIGEST, DAILY_DIGEST, WEEKLY_DIGEST, MONTHLY_DIGEST = range(4)

# Digest options.
DIGEST_MAP = OrderedDict([
    (NO_DIGEST, 'Never'),
    (DAILY_DIGEST, 'Daily'),
    (WEEKLY_DIGEST, 'Weekly'),
    (MONTHLY_DIGEST, 'Monthly')
])

DIGEST_CHOICES = DIGEST_MAP.items()

# Default digest option.
DEFAULT_DIGEST = WEEKLY_DIGEST

# Django allauth configuration
SOCIALACCOUNT_ADAPTER = 'biostar3.middleware.AutoSignupAdapter'

# See the django_allauth docs for what the fields mean.
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