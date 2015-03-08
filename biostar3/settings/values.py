#
# Site specic behaviors.
# These values are loaded first even before the base settings module.
#
# These are located in settings to allow overriding them.
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
POSTS_PER_PAGE = 10

# How many minutes until a post view from an IP is counted again.
POST_VIEW_INTERVAL = 5

# The CSS classes associated with the Django messages framework.
MESSAGE_TAGS = {
    10: 'info', 20: 'info', 25: 'success', 30: 'warning', 40: 'error',
}

# Sort values for userlist.
USER_SORT_BY_VISIT, USER_SORT_BY_REP = "visit", "reputation"
USER_SORT_BY_JOIN, USER_SORT_BY_ACTIVITY = "join", "activity"

# Interface representation for user list sort values.
USER_SORT_CHOICES = [
    (USER_SORT_BY_VISIT, "recent visit"),
    (USER_SORT_BY_REP, "reputation"),
    (USER_SORT_BY_JOIN, "date joined"),
    (USER_SORT_BY_ACTIVITY, "activity level"),
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
    (SORT_BY_UPDATE, "by activity",),
    (SORT_BY_VIEWS, "by views"),
    (SORT_BY_SUBS, "by followers"),
    (SORT_BY_ANSWERS, "by answers"),
    (SORT_BY_BOOKMARKS, "by bookmarks"),
    (SORT_BY_VOTES, "by votes"),
    (SORT_BY_RANK, "by rank"),
    (SORT_BY_CREATION, "by creation"),
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
LOCAL_MESSAGE, EMAIL_MESSAGE, NO_MESSAGES, SMART_MESSAGES, ALL_MESSAGES = range(5)

# The mapping from a messaging type to a readable word.
MESSAGE_CHOICES = [
    (SMART_MESSAGES, "Smart mode"),
    (LOCAL_MESSAGE, "Local messages"),
    (EMAIL_MESSAGE, "Email messages"),
    (ALL_MESSAGES, "Mail List mode"),
    (NO_MESSAGES, "No Messages"),
]

# Default messaging value for a new user.
MESSAGE_DEFAULT = SMART_MESSAGES

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