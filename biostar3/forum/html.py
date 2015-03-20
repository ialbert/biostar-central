"""
Handles the output formatting and conversions
"""
import bleach, re, logging, requests
from django.conf import settings
from markdown2 import markdown
from html5lib.tokenizer import HTMLTokenizer

logger = logging.getLogger('biostar')

# Tags that are allowed by default
ALLOWED_TAGS = bleach.ALLOWED_TAGS + settings.ALLOWED_TAGS
ALLOWED_STYLES = bleach.ALLOWED_STYLES + settings.ALLOWED_STYLES
ALLOWED_ATTRIBUTES = dict(bleach.ALLOWED_ATTRIBUTES)
ALLOWED_ATTRIBUTES.update(settings.ALLOWED_ATTRIBUTES)

# Moderators may use more tags and styles
TRUSTED_TAGS = ALLOWED_TAGS + settings.TRUSTED_TAGS
TRUSTED_STYLES = ALLOWED_STYLES + settings.TRUSTED_STYLES
TRUSTED_ATTRIBUTES = dict(ALLOWED_ATTRIBUTES).update(settings.TRUSTED_ATTRIBUTES)

# Patterns that will be recognized and embedded into the posts as links.
USER_PATTERN = r"http(s)?://%s/u/(?P<uid>(\d+))(/)?" % settings.SITE_DOMAIN
POST_PATTERN1 = r"http(s)?://%s/p/(?P<uid>(\d+))(/)?" % settings.SITE_DOMAIN
POST_PATTERN2 = r"http(s)?://%s/p/\d+/\#(?P<uid>(\d+))(/)?" % settings.SITE_DOMAIN

# Matches gists that may be embeded.
GIST_PATTERN = r"https://gist.github.com/(?P<uid>([\w/]+))"

# Matches Youtube video links.
YOUTUBE_PATTERN = r"http(s)?://www.youtube.com/watch\?v=(?P<uid>(\w+))(/)?"

# Twitter: tweets to embed.
TWITTER_PATTERN = r"http(s)?://twitter.com/\w+/status(es)?/(?P<uid>([\d]+))"


def get_embedded_youtube(uid):
    return '<iframe width="420" height="315" src="//www.youtube.com/embed/%s" frameborder="0" allowfullscreen></iframe>' % uid


def get_embedded_gist(uid):
    return '<script src="https://gist.github.com/%s.js"></script>' % uid


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

# Compile the patterns into regular expressions.
USER_RE = re.compile(USER_PATTERN)
POST_RE1 = re.compile(POST_PATTERN1)
POST_RE2 = re.compile(POST_PATTERN2)
GIST_RE = re.compile(GIST_PATTERN)
YOUTUBE_RE = re.compile(YOUTUBE_PATTERN)
TWITTER_RE = re.compile(TWITTER_PATTERN)


def strip_tags(text):
    "Strip html tags from text"
    result = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
    return result


def clean(text):
    "Sanitize text with no other substitutions"
    result = bleach.clean(text, tags=ALLOWED_TAGS, attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)
    return result


def sanitize(text, user):
    "Sanitize text and expand links to match content"

    if not text.strip():
        # No content there.
        return ""

    # Avoid circular imports.
    from biostar3.forum.models import User, Post

    def internal_links(attrs, new=False):
        """Creates links to internal content"""
        try:

            # Don't resolve the link if a user has already specified a text for it.
            href, _text = attrs['href'], attrs['_text']
            if href != _text:
                return attrs

            # Find and match post URL patterns
            post_patt = POST_RE1.search(href) or POST_RE2.search(href)
            if post_patt:
                uid = post_patt.group("uid")
                attrs['_text'] = Post.objects.get(id=uid).title

            # Find an match user URL patterns.
            user_patt = USER_RE.search(href)
            if user_patt:
                uid = user_patt.group("uid")
                attrs['_text'] = User.objects.get(id=uid).name

        except Exception, exc:
            # This function is a convenience feature.
            # Let's not let it crash the whole post parsing.
            logger.error(exc)

        return attrs

    def require_protocol(attrs, new=False):
        "Linkify only if protocols are present"

        if new:
            href, _text = attrs['href'], attrs['_text']
            if href != _text:
                # This has already been linkified.
                return attrs

            # Don't linkify links with no protocols.
            if href[:4] not in ('http', 'ftp:'):
                return None

        return attrs

    # The functions that will be applied when linkifying
    callbacks = [internal_links, require_protocol]

    # Staff may use more dangerous HTML objects.
    if user.is_staff:
        tags, attrs, styles = TRUSTED_TAGS, TRUSTED_ATTRIBUTES, TRUSTED_STYLES
    else:
        tags, attrs, styles = ALLOWED_TAGS, ALLOWED_ATTRIBUTES, ALLOWED_STYLES

    # Sanitize html input.
    html = bleach.clean(text, tags=tags, attributes=attrs, styles=styles)

    try:
        # Apply the markdown transformation.
        # We'll protect against library crashes by a generic Exception catch.
        html = markdown(html, extras=["fenced-code-blocks", "code-friendly", "nofollow", "spoiler"])
    except Exception, exc:
        logger.error('crash during markdown conversion: %s' % exc)
        html = html

    # Find embeddable patterns.
    html = embed_links(html)

    try:
        # Turn links into urls. We had bleach.linkify crash very rarely so we'll trap that.
        # We use a more lenient tokenizer since the content is already cleaned.
        html = bleach.linkify(html, callbacks=callbacks, skip_pre=True, tokenizer=HTMLTokenizer)
    except Exception, exc:
        logger.error('crash during bleach linkify: %s' % exc)
        html = html

    # Strip whitespace.
    html = html.strip()

    return html


def embed_links(text):
    targets = [
        (GIST_RE, get_embedded_gist),
        (YOUTUBE_RE, get_embedded_youtube),
        (TWITTER_RE, get_embedded_tweet),
    ]

    for regex, func in targets:
        patt = regex.search(text)
        if patt:
            uid = patt.group("uid")
            new_value = func(uid)
            text = regex.sub(new_value, text)

    return text

