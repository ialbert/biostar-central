import re
import bleach
import logging
import uuid
import  requests
from datetime import datetime, timedelta


from django.template import loader, Context
from django.utils.timezone import utc

# Allowed html content.
ALLOWED_TAGS = "p div br code pre h1 h2 h3 h4 hr span s sub sup b i img strong strike em underline super table thead tr th td tbody".split()
ALLOWED_STYLES = 'color font-weight background-color width height'.split()
DEFAULT_ATTRIBUTES = {
    '*': ['class', 'style'],
    'a': ['href', 'rel'],
    'img': ['src', 'alt', 'width', 'height'],
    'table': ['border', 'cellpadding', 'cellspacing'],

}

ALLOWED_TAGS = bleach.ALLOWED_TAGS + ALLOWED_TAGS
ALLOWED_STYLES = bleach.ALLOWED_STYLES + ALLOWED_STYLES
ALLOWED_ATTRIBUTES = dict(bleach.ALLOWED_ATTRIBUTES)
ALLOWED_ATTRIBUTES.update(DEFAULT_ATTRIBUTES)


def user_pattern(base_site):
    user = fr"^http(s)?://{base_site}/accounts/profile/(?P<uid>(\w+))(/)?$"
    return user


def post_pattern(base_site):
    post_pattern = fr"^http(s)?://{base_site}/p/(?P<uid>(\w+))(/)?$"
    return post_pattern


# Matches gists that may be embeded.
GIST_PATTERN = r"^https://gist.github.com/(?P<uid>([\w/]+))"
# Matches Youtube video links.
YOUTUBE_PATTERN1 = r"^http(s)?://www.youtube.com/watch\?v=(?P<uid>([\w-]+))(/)?"
YOUTUBE_PATTERN2 = r"https://www.youtube.com/embed/(?P<uid>([\w-]+))(/)?"
YOUTUBE_PATTERN3 = r"https://youtu.be/(?P<uid>([\w-]+))(/)?"

# Twitter: tweets to embed.
TWITTER_PATTERN = r"http(s)?://twitter.com/\w+/status(es)?/(?P<uid>([\d]+))"

logger = logging.getLogger(__name__)


def clean(text):
    """
    Sanitize text with no other substitutions
    """
    html = bleach.clean(text, tags=ALLOWED_TAGS,
                        attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)
    return html


def get_embedded_tweet(tweet_id):
    """
    Get the HTML code with the embedded tweet.
    It requires an API call at https://api.twitter.com/1/statuses/oembed.json as documented here:
    https://dev.twitter.com/docs/embedded-tweets - section "Embedded Tweets for Developers"
    https://dev.twitter.com/docs/api/1/get/statuses/oembed

    Params:
    tweet_id -- a tweet's numeric id like 2311234267 for the tweet at
    https://twitter.com/Linux/status/2311234267
    """
    try:
        response = requests.get("https://api.twitter.com/1/statuses/oembed.json?id={}".format(
            tweet_id))
        return response.json()['html']
    except:
        return ''


def get_target_text(pattern):



    return pattern_map.get(pattern)


def fixcase(name):
    return name.upper() if len(name) == 1 else name.lower()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def split_tags(text):

    capitalize = lambda txt: txt.upper() if len(txt) == 1 else txt
    return [capitalize(x) for x in text.split(",") if len(x)]


def render(name, **kwds):
    "Helper function to render a template"
    tmpl = loader.get_template(name)
    page = tmpl.render(kwds)
    return page


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def strip_tags(text):
    "Strip html tags from text"
    text = bleach.clean(text, tags=[], attributes={}, styles=[], strip=True)
    return text


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


def time_ago(date):

    if not date:
        return ''
    delta = now() - date
    if delta < timedelta(minutes=1):
        return 'just now'
    elif delta < timedelta(hours=1):
        unit = pluralize(delta.seconds // 60, "minute")
    elif delta < timedelta(days=1):
        unit = pluralize(delta.seconds // 3600, "hour")
    elif delta < timedelta(days=30):
        unit = pluralize(delta.days, "day")
    elif delta < timedelta(days=90):
        unit = pluralize(int(delta.days / 7), "week")
    elif delta < timedelta(days=730):
        unit = pluralize(int(delta.days / 30), "month")
    else:
        diff = delta.days / 365.0
        unit = '%0.1f years' % diff
    return "%s ago" % unit