import re
import bleach
import logging
import requests
import mistune
import uuid
from datetime import datetime
#from html5lib.tokenizer import HTMLTokenizer


from . import const
from django.template import loader, Context
from django.utils.timezone import utc

logger = logging.getLogger(__name__)

ALLOWED_TAGS = bleach.ALLOWED_TAGS + const.ALLOWED_TAGS
ALLOWED_STYLES = bleach.ALLOWED_STYLES + const.ALLOWED_STYLES
ALLOWED_ATTRIBUTES = dict(bleach.ALLOWED_ATTRIBUTES)
ALLOWED_ATTRIBUTES.update(const.ALLOWED_ATTRIBUTES)

# Matching patterns will be filled in with post title or user name
#TODO: Change this to match
USER_PATTERN = r"^http(s)?://%s/u/(?P<uid>(\d+))(/)?$" % const.SITE_DOMAIN
POST_PATTERN1 = r"^http(s)?://%s/p/(?P<uid>(\d+))(/)?$" % const.SITE_DOMAIN
POST_PATTERN2 = r"^http(s)?://%s/p/\d+/\#(?P<uid>(\d+))(/)?$" % const.SITE_DOMAIN

# Matches gists that may be embeded.
GIST_PATTERN = r"^https://gist.github.com/(?P<uid>([\w/]+))"

# Matches Youtube video links.
YOUTUBE_PATTERN1 = r"^http(s)?://www.youtube.com/watch\?v=(?P<uid>([\w-]+))(/)?"
YOUTUBE_PATTERN2 = r"https://www.youtube.com/embed/(?P<uid>([\w-]+))(/)?"
YOUTUBE_PATTERN3 = r"https://youtu.be/(?P<uid>([\w-]+))(/)?"

# Twitter: tweets to embed.
TWITTER_PATTERN = r"http(s)?://twitter.com/\w+/status(es)?/(?P<uid>([\d]+))"

USER_RE = re.compile(USER_PATTERN)
POST_RE1 = re.compile(POST_PATTERN1)
POST_RE2 = re.compile(POST_PATTERN2)
GIST_RE = re.compile(GIST_PATTERN)
YOUTUBE_RE1 = re.compile(YOUTUBE_PATTERN1)
YOUTUBE_RE2 = re.compile(YOUTUBE_PATTERN2)
YOUTUBE_RE3 = re.compile(YOUTUBE_PATTERN3)

TWITTER_RE = re.compile(TWITTER_PATTERN)


def fixcase(name):
    return name.upper() if len(name) == 1 else name.lower()


def now():
    return datetime.utcnow().replace(tzinfo=utc)

def split_tags(text):

    capitalize = lambda txt: txt.upper() if len(txt)==1 else txt

    return [capitalize(x) for x in text.split(",") if len(x)]

def clean(text):
    "Sanitize text with no other substitutions"
    html = bleach.clean(text, tags=ALLOWED_TAGS,
                        attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)
    return html

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


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


def strip_tags(text):
    "Strip html tags from text"
    text = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
    return text


def render(name, **kwds):
    "Helper function to render a template"
    tmpl = loader.get_template(name)
    cont = Context(kwds)
    page = tmpl.render(cont)
    return page


def test():
    pass


if __name__ == '__main__':
    test()